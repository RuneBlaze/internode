use clap::ArgEnum;
use rayon::string;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum Mode {
    Support,
    Internode,
    NLength,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum ImputeMethod {
    Upgma,
    BalMENNI,
    BalMESPR,
}

pub struct UstarConfig {
    pub upper_bound: f64,
    pub lower_bound: f64,
    pub mode: Mode,
}

impl Default for UstarConfig {
    fn default() -> Self {
        UstarConfig {
            upper_bound: 1.0,
            lower_bound: 0.0,
            mode: Mode::Support,
        }
    }
}

#[derive(Debug)]
pub struct TreeCollection {
    pub taxon_set: TaxonSet,
    pub trees: Vec<Tree>,
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

#[derive(Debug)]
pub struct TaxonSet {
    pub to_id: HashMap<String, usize>,
    pub names: Vec<String>,
    last: usize,
}

impl TaxonSet {
    pub fn request(&mut self, taxon_name: String) -> usize {
        *self.to_id.entry(taxon_name.clone()).or_insert_with(|| {
            self.names.push(taxon_name);
            self.last += 1;
            self.last - 1
        })
    }

    pub fn retrieve(&self, taxon_name: &str) -> usize {
        *self.to_id.get(taxon_name).unwrap()
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
            Ok(TreeCollection { taxon_set, trees })
        } else {
            Err("Could not read file")
        }
    }

    pub fn ngenes(&self) -> usize {
        self.trees.len()
    }

    pub fn ntaxa(&self) -> usize {
        self.taxon_set.len()
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
    pub root: usize,
}

impl Tree {
    pub fn rooted_star(size: usize) -> Self {
        let mut taxa = vec![-1; 2 * size - 1];
        for i in 0..size {
            taxa[i] = i as i32;
        }
        let parents = vec![-1; 2 * size - 1];
        let support = vec![0.0; 2 * size - 1];
        let lengths = vec![0.0; 2 * size - 1];
        let firstchild = vec![-1; 2 * size - 1];
        let nextsib = vec![-1; 2 * size - 1];
        let childcount = vec![0; 2 * size - 1];
        let fake_root = true;
        Tree {
            taxa,
            parents,
            support,
            lengths,
            firstchild,
            nextsib,
            childcount,
            fake_root,
            root: 0,
        }
    }

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
        node == self.root
    }

    pub fn topology_newick(&self, taxon_set: &TaxonSet) -> String {
        let mut string_rep: Vec<String> = vec![String::new(); self.taxa.len()];
        for node in self.postorder() {
            if self.is_leaf(node) {
                string_rep[node] = taxon_set.names[self.taxa[node] as usize].clone();
            } else {
                let mut out = String::new();
                out.push_str("(");
                for c in self.children(node) {
                    out.push_str(&string_rep[c]);
                    out.push_str(",");
                }
                out.pop();
                out.push_str(")");
                string_rep[node] = out;
            }
        }
        string_rep[self.root].push_str(";");
        string_rep.swap_remove(self.root)
    }
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
            if !ls.is_empty() {
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
                let rg = config.upper_bound - config.lower_bound;
                support[n] = ((s - config.lower_bound) / rg).max(0.0);
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

    Tree {
        taxa,
        parents,
        support,
        lengths,
        firstchild,
        nextsib,
        childcount,
        fake_root,
        root: 0,
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
            tree,
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
        s1.push(tree.root);
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
