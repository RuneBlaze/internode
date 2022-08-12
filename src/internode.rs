#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]

use crate::tree::*;
use ndarray::prelude::*;
use ndarray::Array;
use rayon::prelude::*;
use std::cell::RefCell;
use std::ffi::CString;
use std::mem::size_of;
use std::ptr;
use std::sync::Arc;
use thread_local::ThreadLocal;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[derive(Debug)]
pub struct UstarState {
    pub dm: Array<f64, Ix2>,
    pub mask: Array<u32, Ix2>,
    pub dim: usize,
    pub minted: bool,
    pub temp: Option<Array<f64, Ix2>>,
    pub norm_factor: f64,
    pub has_missing: bool,
}

impl UstarState {
    pub fn from_taxon_set(taxon_set: &TaxonSet, config: &UstarConfig) -> Self {
        let n = taxon_set.len();
        let dm = Array::<f64, _>::zeros((n, n).f());
        let mask = Array::<u32, _>::zeros((n, n).f());
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
            has_missing: false,
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
        state
    }

    pub fn from_tree_collection_par(
        tree_collection: &TreeCollection,
        config: &UstarConfig,
        nthreads: usize,
    ) -> Self {
        let tls = Arc::new(ThreadLocal::new());
        // FIXME: this is really a hack for efficient parallelization
        let chunk_size_bound = match tree_collection.ntaxa() {
            n if n <= 50 => 10000usize,
            n if n <= 500 => 2000usize,
            n if n <= 1000 => 500usize,
            _ => 200usize,
        };
        let chunk_size = (tree_collection.trees.len() / nthreads + 1).max(chunk_size_bound);
        let _ = &tree_collection
            .trees
            .par_chunks(chunk_size)
            .for_each(|trees| {
                let tls2 = tls.clone();
                let state = tls2.get_or(|| {
                    RefCell::new(UstarState::from_taxon_set(
                        &tree_collection.taxon_set,
                        config,
                    ))
                });
                let mut borrowed = state.borrow_mut();
                for t in trees {
                    add_to_matrix(&mut borrowed, t, config.mode);
                }
            });
        let mut state = UstarState::from_taxon_set(&tree_collection.taxon_set, config);
        Arc::try_unwrap(tls).unwrap().into_iter().for_each(|s| {
            state.add_from(&s.borrow());
        });
        state
    }

    pub fn flatten(&mut self) {
        for i in 0..self.dim {
            for j in (i + 1)..self.dim {
                if self.mask[[i, j]] <= 0 {
                    self.has_missing = true;
                } else {
                    self.dm[[i, j]] /= self.mask[[i, j]] as f64;
                }
            }
        }
        self.minted = true;
    }

    pub fn raw_tree(&mut self, taxon_set: &TaxonSet) -> String {
        run_fastme(taxon_set, &self.dm, &FastMEConfig::default())
    }

    pub fn to_tree(&mut self, taxon_set: &TaxonSet) -> String {
        if !self.minted {
            self.flatten();
        }
        self.raw_tree(taxon_set)
    }

    pub fn add_from(&mut self, rhs: &UstarState) {
        rhs.dm.indexed_iter().for_each(|((i, j), v)| {
            self.dm[[i, j]] += v;
        });
        rhs.mask.indexed_iter().for_each(|((i, j), v)| {
            self.mask[[i, j]] += v;
        });
    }
}

pub fn add_to_matrix(state: &mut UstarState, tree: &Tree, mode: Mode) {
    // a straightforward translation of the Treeswift logic
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
                            state.mask[[l, r]] += 1;
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

// FIXME: use trait to DRY
pub fn impute_matrix(state: &mut UstarState, tree: &Tree, mode: Mode) {
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
                            if state.mask[[l, r]] <= 0 {
                                state.dm[[l, r]] = dist;
                            }
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
    let mut max_dis: f64 = 0.0;
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
                            state.mask[[l, r]] += 1;
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

pub struct FastMEConfig {
    use_nni: bool,
    use_spr: bool,
}

impl FastMEConfig {
    pub fn new(use_nni: bool, use_spr: bool) -> FastMEConfig {
        FastMEConfig { use_nni, use_spr }
    }
}

impl Default for FastMEConfig {
    fn default() -> Self {
        FastMEConfig {
            use_nni: true,
            use_spr: true,
        }
    }
}

pub fn run_fastme(
    taxon_set: &TaxonSet,
    dm: &Array<f64, Ix2>,
    fastme_config: &FastMEConfig,
) -> String {
    let size = taxon_set.len() as i32;
    unsafe {
        let mut A = initDoubleMatrix(2 * size - 2);
        let D = initDoubleMatrix(2 * size - 2);
        fillZeroMatrix(&mut A as *mut *mut *mut f64, 2 * size - 2);
        for i in 0..size {
            let ptr = D.offset(i as isize);
            *ptr = mCalloc(size, size_of::<f64>() as u64) as *mut f64;
        }
        for i in 0..size {
            for j in (i + 1)..size {
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
        options.use_NNI = if fastme_config.use_nni { 1 } else { 0 };
        options.use_SPR = if fastme_config.use_spr { 1 } else { 0 };
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
        let translated = translate_newick(taxon_set, &s);
        freeMatrix(A, 2 * size - 2);
        freeMatrix(D, 2 * size - 2); // there are probably other leaks present
        translated
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
            if !ls.is_empty() {
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
            if !buf.ends_with(')') {
                let leaf_ix = ts.parse::<usize>().unwrap();
                buf.push_str(&taxon_set.names[leaf_ix]);
            } else {
                buf.push_str(&ts);
            }
        }
    }
    buf
}

#[cfg(test)]
pub mod tests {
    use std::path::PathBuf;
    pub fn avian_tree() -> PathBuf {
        let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        d.push("resources/test");
        d.push("avian.tre");
        d
    }
}
