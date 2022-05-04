#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]

use clap::ArgEnum;
use ndarray::AssignElem;
use ndarray::prelude::*;
use ndarray::Array;
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::collections::VecDeque;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufRead};
use std::mem::size_of;
use std::ops::AddAssign;
use std::path::Path;
use std::ptr;
use std::sync::Arc;
use std::{borrow::Borrow, cmp::max, collections::HashMap};
use rayon::prelude::*;
use thread_local::ThreadLocal;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[derive(Debug)]
pub struct TaxonSet {
    pub to_id: HashMap<String, usize>,
    pub names: Vec<String>,
    last: usize,
}
// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

impl TaxonSet {
    pub fn request(&mut self, taxon_name: String) -> usize {
        self.to_id
            .entry(taxon_name.clone())
            .or_insert_with(|| {
                self.names.push(taxon_name);
                self.last += 1;
                self.last - 1
            })
            .clone()
    }

    pub fn retreive(&self, taxon_name: String) -> usize {
        self.to_id.get(&taxon_name).unwrap().clone()
    }

    pub fn new() -> Self {
        TaxonSet {
            to_id: HashMap::new(),
            names: Vec::new(),
            last: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.last
    }
}

#[derive(Debug)]
pub struct UstarState {
    pub dm: Array<f64, Ix2>,
    pub mask: Array<f64, Ix2>,
    pub dim: usize,
    pub minted: bool,
    pub temp: Option<Array<f64, Ix2>>,
    pub norm_factor : f64,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum Mode {
    Support,
    Internode,
    NLength,
}

pub struct UstarConfig {
    pub max_support: f64,
    pub normalizer: f64,
    pub mode: Mode,
}

impl Default for UstarConfig {
    fn default() -> Self {
        UstarConfig {
            max_support: 1.0,
            normalizer: 0.0,
            mode: Mode::Support,
        }
    }
}

impl UstarState {
    pub fn from_taxon_set(taxon_set: &TaxonSet, config: &UstarConfig) -> Self {
        let n = taxon_set.len();
        let dm = Array::<f64, _>::zeros((n, n).f());
        let mask = Array::<f64, _>::zeros((n, n).f());
        UstarState {
            dm,
            mask,
            dim: n,
            minted: false,
            temp: match config.mode {
                Mode::NLength => Some(Array::<f64, _>::zeros((n, n).f())),
                _ => None,
            },
            norm_factor: -1.0,
        }
    }

    pub fn from_tree_collection(tree_collection: &TreeCollection, config: &UstarConfig) -> Self {
        let mut state = UstarState::from_taxon_set(&tree_collection.taxon_set, config);
        if config.mode == Mode::NLength {
            state.temp = Some(Array::<f64, _>::zeros((state.dim, state.dim).f()));
            for t in &tree_collection.trees {
                add_to_matrix_with_temp(&mut state, t, config.mode);
            }
        } else {
            for t in &tree_collection.trees {
                add_to_matrix(&mut state, t, config.mode);
            }
        }
        return state;
    }

    pub fn from_tree_collection_par(tree_collection: &TreeCollection, config: &UstarConfig) -> Self {
        let tls = Arc::new(ThreadLocal::new());
        let _ = &tree_collection.trees.par_iter().for_each(|t| {
            let state = tls.get_or(|| 
                RefCell::new(UstarState::from_taxon_set(&tree_collection.taxon_set, config)));
            add_to_matrix(&mut state.borrow_mut(), t, config.mode);
        });
        let mut state = UstarState::from_taxon_set(&tree_collection.taxon_set, config);
        Arc::try_unwrap(tls).unwrap().into_iter().for_each(|s| {
            state.add_from(&s.borrow());
        });
        return state;
    }

    pub fn flatten(&mut self) {
        for i in 0..self.dim {
            for j in (i + 1)..self.dim {
                self.dm[[i, j]] /= self.mask[[i, j]];
            }
        }
        self.minted = true;
    }

    pub fn raw_tree(&mut self, taxon_set: &TaxonSet) -> String {
        run_fastme(taxon_set, &self.dm)
    }

    pub fn to_tree(&mut self, taxon_set: &TaxonSet) -> String {
        if !self.minted {
            self.flatten();
        }
        self.raw_tree(taxon_set)
    }

    pub fn add_from(&mut self, rhs : &UstarState) {
        rhs.dm.indexed_iter().for_each(|((i, j), v)| {
            self.dm[[i, j]] += v;
        });
        rhs.mask.indexed_iter().for_each(|((i, j), v)| {
            self.mask[[i, j]] += v;
        });
    }
}

pub fn add_to_matrix(state: &mut UstarState, tree: &Tree, mode: Mode) {
    // a straightforward translation of the treeswift logic
    // sparse vector of distances
    let mut leaf_dists = Vec::<Vec<(usize, f64)>>::new();
    leaf_dists.resize(tree.taxa.len(), Vec::new());
    for node in tree.postorder() {
        if tree.is_leaf(node) {
            leaf_dists[node].push((node, 0.0));
        } else {
            let mut calculated_root = false;
            for c in tree.children(node) {
                if tree.is_root(node) && tree.fake_root {
                    if calculated_root {
                        continue;
                    }
                    calculated_root = true;
                }
                for e in leaf_dists.get_mut(c).unwrap() {
                    e.1 += match mode {
                        Mode::Support => tree.support[c],
                        Mode::Internode => 1.0,
                        Mode::NLength => tree.lengths[c],
                    };
                }
            }

            let node_children: Vec<usize> = tree.children(node).collect();
            for c1 in 0..(node_children.len() - 1) {
                let leaves_c1 = leaf_dists.get(node_children[c1]).unwrap();
                for c2 in (c1 + 1)..(node_children.len()) {
                    let leaves_c2 = leaf_dists.get(node_children[c2]).unwrap();
                    for i in 0..(leaves_c1.len()) {
                        for j in 0..(leaves_c2.len()) {
                            let (u, ud) = leaves_c1[i];
                            let (v, vd) = leaves_c2[j];
                            let dist = ud + vd;
                            let u_leaf = tree.taxa[u] as usize;
                            let v_leaf = tree.taxa[v] as usize;
                            let l = std::cmp::min(u_leaf, v_leaf);
                            let r = std::cmp::max(u_leaf, v_leaf);
                            state.dm[[l, r]] += dist;
                            state.mask[[l, r]] += 1.0;
                        }
                    }
                }
            }

            for (i, e) in tree.children(node).enumerate() {
                if i == 0 {
                    leaf_dists.swap(node, tree.firstchild[node] as usize);
                } else {
                    leaf_dists.push(vec![]);
                    let mut v = leaf_dists.swap_remove(e);
                    leaf_dists[node].append(&mut v);
                }
            }
        }
    }
}

// FIXME: this is duplicating code
// used only when mode is NLength
pub fn add_to_matrix_with_temp(state: &mut UstarState, tree: &Tree, _: Mode) {
    let temp = state.temp.as_mut().unwrap();
    temp.fill(0.0);
    let mut leaf_dists = Vec::<Vec<(usize, f64)>>::new();
    leaf_dists.resize(tree.taxa.len(), Vec::new());
    let mut max_dis : f64 = 0.0;
    for node in tree.postorder() {
        if tree.is_leaf(node) {
            leaf_dists[node].push((node, 0.0));
        } else {
            let mut calculated_root = false;
            for c in tree.children(node) {
                if tree.is_root(node) && tree.fake_root {
                    if calculated_root {
                        continue;
                    }
                    calculated_root = true;
                }
                for e in leaf_dists.get_mut(c).unwrap() {
                    e.1 += tree.lengths[c];
                }
            }

            let node_children: Vec<usize> = tree.children(node).collect();
            for c1 in 0..(node_children.len() - 1) {
                let leaves_c1 = leaf_dists.get(node_children[c1]).unwrap();
                for c2 in (c1 + 1)..(node_children.len()) {
                    let leaves_c2 = leaf_dists.get(node_children[c2]).unwrap();
                    for i in 0..(leaves_c1.len()) {
                        for j in 0..(leaves_c2.len()) {
                            let (u, ud) = leaves_c1[i];
                            let (v, vd) = leaves_c2[j];
                            let dist = ud + vd;
                            let u_leaf = tree.taxa[u] as usize;
                            let v_leaf = tree.taxa[v] as usize;
                            let l = std::cmp::min(u_leaf, v_leaf);
                            let r = std::cmp::max(u_leaf, v_leaf);
                            temp[[l, r]] += dist;
                            max_dis = max_dis.max(dist);
                            state.mask[[l, r]] += 1.0;
                        }
                    }
                }
            }

            for (i, e) in tree.children(node).enumerate() {
                if i == 0 {
                    leaf_dists.swap(node, tree.firstchild[node] as usize);
                } else {
                    leaf_dists.push(vec![]);
                    let mut v = leaf_dists.swap_remove(e);
                    leaf_dists[node].append(&mut v);
                }
            }
        }
    }
    if state.norm_factor <= 0.0 {
        state.norm_factor = max_dis;
    }
    if max_dis <= 0.0 {
        return;
    }
    max_dis /= state.norm_factor;
    temp.indexed_iter().for_each(|((i, j), v)| {
        state.dm[[i, j]] += v / max_dis;
    });
}

#[derive(Debug)]
pub struct TreeCollection {
    pub taxon_set: TaxonSet,
    pub trees: Vec<Tree>,
}

impl TreeCollection {
    pub fn new() -> Self {
        TreeCollection {
            taxon_set: TaxonSet::new(),
            trees: Vec::new(),
        }
    }

    pub fn from_newick<P>(filename: P, config: &UstarConfig) -> Result<Self, &'static str>
    where
        P: AsRef<Path>,
    {
        let mut trees: Vec<Tree> = vec![];
        let mut taxon_set = TaxonSet::new();
        if let Ok(lines) = read_lines(filename) {
            for line in lines {
                if let Ok(newick) = line {
                    let parsed = parse_newick(&mut taxon_set, newick.as_str(), config);
                    trees.push(parsed);
                } else {
                    return Err("Error reading file");
                }
            }
            return Ok(TreeCollection { taxon_set, trees });
        } else {
            return Err("Could not read file");
        }
    }
}

#[derive(Debug)]
pub struct Tree {
    pub taxa: Vec<i32>,
    pub parents: Vec<i32>,
    pub support: Vec<f64>, // branch support
    pub lengths: Vec<f64>, // branch lengths
    pub firstchild: Vec<i32>,
    pub nextsib: Vec<i32>,
    pub childcount: Vec<u32>,
    pub fake_root: bool,
}

impl Tree {
    pub fn children(&self, node: usize) -> ChildrenIterator {
        ChildrenIterator::new(self, node)
    }

    pub fn postorder(&self) -> PostorderIterator {
        PostorderIterator::new(self)
    }

    pub fn is_leaf(&self, node: usize) -> bool {
        if self.childcount[node] == 0 {
            return true;
        }
        false
    }

    pub fn is_root(&self, node: usize) -> bool {
        return node == 0;
    }
}

pub struct PostorderIterator {
    // s1 : Vec<usize>,
    s2: Vec<usize>,
}

pub struct ChildrenIterator<'a> {
    tree: &'a Tree,
    current: i32,
}

impl<'a> ChildrenIterator<'a> {
    pub fn new(tree: &'a Tree, node: usize) -> Self {
        ChildrenIterator {
            tree: tree,
            current: tree.firstchild[node],
        }
    }
}

impl<'a> Iterator for ChildrenIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == -1 {
            None
        } else {
            let res = self.current as usize;
            self.current = self.tree.nextsib[self.current as usize];
            Some(res)
        }
    }
}

impl PostorderIterator {
    pub fn new(tree: &Tree) -> Self {
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        s1.push(0usize);
        while let Some(n) = s1.pop() {
            s2.push(n);
            tree.children(n).for_each(|c| s1.push(c));
        }
        PostorderIterator {
            // s1,
            s2,
        }
    }
}

impl Iterator for PostorderIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.s2.pop()
    }
}

pub fn run_fastme(taxon_set: &TaxonSet, dm: &Array<f64, Ix2>) -> String {
    let size = taxon_set.len() as i32;
    unsafe {
        let mut A = initDoubleMatrix(2 * size - 2);
        let mut D = initDoubleMatrix(2 * size - 2);
        fillZeroMatrix(&mut A as *mut *mut *mut f64, 2 * size - 2);
        for i in 0..size {
            let ptr = D.offset(i as isize);
            *ptr = mCalloc(size, size_of::<f64>() as u64) as *mut f64;
        }
        for i in 0..size {
            for j in (i + 1)..size {
                // let row = D.offset(i as isize);
                // let r2 = *row;
                // let ptr = r2.offset(j as isize);

                let r1 = (*(D.offset(i as isize))).offset(j as isize);
                let r2 = (*(D.offset(j as isize))).offset(i as isize);
                let l = i as usize;
                let r = j as usize;
                *r1 = dm[[l, r]];
                *r2 = dm[[l, r]];
            }
        }

        let mut options = Options::default();
        Set_Defaults_Input(&mut options);
        options.use_SPR = 1;
        options.use_NNI = 1;
        options.method = TaxAddBAL as i32;
        options.NNI = BALNNI as i32;

        let mut species = set::default(); // terrible naming! Bad FastME!
        species.firstNode = ptr::null_mut();
        species.secondNode = ptr::null_mut();

        for (i, _) in taxon_set.names.iter().enumerate() {
            let name = CString::new(i.to_string()).unwrap();
            let mut v = makeNode(name.as_ptr(), -1);
            (*v).index2 = i as i32;
            addToSet(v, &mut species as *mut set);
        }

        let mut t = ComputeTree(&mut options, D, A, &mut species, size, 8);
        let mut nniCount: i32 = 0;
        let mut sprCount: i32 = 0;
        t = ImproveTree(
            &mut options,
            t,
            D,
            A,
            &mut nniCount,
            &mut sprCount,
            options.fpO_stat_file,
        );
        let mut tree_output = vec![0u8; (size << 10) as usize]; //Vec::<u8>::with_capacity((size << 10) as usize);
        NewickPrintTreeStr(t, tree_output.as_mut_ptr() as *mut i8, 2);
        let result = CString::from_vec_unchecked(tree_output);
        let s = result.to_str().unwrap().to_owned();
        return translate_newick(taxon_set, &s);
    }
}

pub fn translate_newick(taxon_set: &TaxonSet, newick: &str) -> String {
    let mut buf = String::new();
    let mut chars = newick.chars().fuse().peekable();
    while let Some(c) = chars.next() {
        if c == ';' {
            buf.push(';');
            break;
        } else if c == '(' {
            buf.push('(');
        } else if c == ')' {
            buf.push(')');
        } else if c == ',' {
            buf.push(',');
        } else if c == ':' {
            let mut ls = "".to_string();
            loop {
                match chars.peek() {
                    Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ls.push(chars.next().unwrap());
                    }
                }
            }
            if ls.len() > 0 {
                buf.push(':');
                buf.push_str(&ls);
            }
        } else {
            let mut ts = c.to_string();
            loop {
                match chars.peek() {
                    Some(':') | Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ts.push(chars.next().unwrap());
                    }
                }
            }
            if !buf.ends_with(")") {
                let leaf_ix = ts.parse::<usize>().unwrap();
                buf.push_str(&taxon_set.names[leaf_ix]);
            } else {
                buf.push_str(&ts);
            }
        }
    }
    return buf;
}

pub fn parse_newick(taxon_set: &mut TaxonSet, newick: &str, config: &UstarConfig) -> Tree {
    let mut taxa: Vec<i32> = vec![-42];
    let mut parents: Vec<i32> = vec![0];
    let mut support: Vec<f64> = vec![-1.0];
    let mut lengths: Vec<f64> = vec![-1.0];
    let mut childcount: Vec<u32> = vec![0];
    let mut firstchild: Vec<i32> = vec![-1];
    let mut nextsib: Vec<i32> = vec![-1];
    // we just reuse TreeSwift's logic
    let mut n: usize = 0; // the current node
    let mut chars = newick.chars().fuse().peekable();
    while let Some(c) = chars.next() {
        if c == ';' {
            break;
        } else if c == '(' {
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            support.push(0.0);
            lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            firstchild[n] = (taxa.len() - 1) as i32;
            n = taxa.len() - 1;
        } else if c == ')' {
            n = parents[n] as usize;
        } else if c == ',' {
            nextsib[n] = (taxa.len()) as i32;
            n = parents[n] as usize;
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            support.push(0.0);
            lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            n = taxa.len() - 1;
        } else if c == ':' {
            let mut ls = "".to_string();
            loop {
                match chars.peek() {
                    Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ls.push(chars.next().unwrap());
                    }
                }
            }
            if ls.len() > 0 {
                lengths[n as usize] = ls.parse::<f64>().unwrap();
            }
        } else {
            let mut ts = c.to_string();
            loop {
                match chars.peek() {
                    Some(':') | Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ts.push(chars.next().unwrap());
                    }
                }
            }
            if childcount[n] == 0 {
                taxa[n] = taxon_set.request(ts) as i32;
                support[n] = 1.0;
            } else {
                let s = ts.parse::<f64>().unwrap();
                let normalizer = config.normalizer;
                support[n] = if normalizer <= 0.0 {
                    s / config.max_support
                } else {
                    (s / config.max_support - config.normalizer) / (1.0 - config.normalizer)
                };
            }
        }
    }

    let mut fake_root = false;
    if childcount[0] == 2 {
        let c = firstchild[0] as usize;
        let c2 = nextsib[c] as usize;
        // then the root is "fake". We need to correct the support values
        fake_root = true;
        // the philosophy here is that this fake edge should only be traversed once
        let supp = support[c].max(support[c2]);
        support[c] = supp;
        support[c2] = supp;

        let length = lengths[c] + lengths[c2];
        lengths[c] = length;
        lengths[c2] = length;
    }
    let tree = Tree {
        taxa,
        parents,
        support,
        lengths,
        firstchild,
        nextsib,
        childcount,
        fake_root,
    };
    return tree;
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::path::PathBuf;

    pub fn avian_tree() -> PathBuf {
        let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        d.push("resources/test");
        d.push("avian.tre");
        d
    }

    #[test]
    fn test_read_tree_collection() {
        let trees = TreeCollection::from_newick(avian_tree(), &UstarConfig::default()).unwrap();
        assert_eq!(48, trees.taxon_set.len());
        assert_eq!(100, trees.trees.len());
    }

    #[test]
    fn test_ustar_state_from_collection() {
        let trees = TreeCollection::from_newick(avian_tree(), &UstarConfig::default()).unwrap();
        let ustar = UstarState::from_tree_collection(&trees, &UstarConfig::default());
        assert!(!ustar.mask.is_empty());
    }
}