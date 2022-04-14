#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]

use std::{collections::HashMap, borrow::Borrow, cmp::max};
use std::mem::{size_of};
use std::ptr;
use std::ffi::CString;
use ndarray::Array;
use ndarray::prelude::*;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[derive(Debug)]
pub struct TaxonSet {
    pub to_id : HashMap<String, usize>,
    pub names : Vec<String>,
    last : usize,
}

impl TaxonSet {
    pub fn request(&mut self, taxon_name : String) -> usize {
        self.to_id.entry(taxon_name.clone()).or_insert_with(|| {
            self.names.push(taxon_name);
            self.last += 1;
            self.last - 1
        }).clone()
    }

    pub fn new() -> Self {
        TaxonSet {
            to_id : HashMap::new(),
            names : Vec::new(),
            last : 0,
        }
    }

    pub fn len(&self) -> usize {
        self.last
    }
}

#[derive(Debug)]
pub struct TreeCollection {
    pub taxon_set : TaxonSet,
    pub trees : Vec<Tree>,
}

impl TreeCollection {
    pub fn new() -> Self {
        TreeCollection {
            taxon_set : TaxonSet::new(),
            trees : Vec::new(),
        }
    }

    pub fn from_newick() -> Self {
        unimplemented!()
    }
}

#[derive(Debug)]
pub struct Tree {
    pub taxa : Vec<i32>,
    pub parents : Vec<i32>,
    pub support : Vec<f64>, // branch support
    pub lengths : Vec<f64>, // branch lengths
    pub firstchild : Vec<i32>,
    pub nextsib : Vec<i32>,
    pub childcount : Vec<u32>,
    pub fake_root : bool,
}

pub fn run_fastme(taxon_set : &TaxonSet, dm : &Array<f64, Ix2>) -> String {
    let size = taxon_set.len() as i32;
    unsafe {
        // println!("size {}", size);
        // println!("initializing double matrix");
        let mut A = initDoubleMatrix(2*size-2);
        let mut D = initDoubleMatrix(2*size-2);
        fillZeroMatrix(&mut A as *mut *mut *mut f64, 2 * size-2);
        for i in 0..size {
            let ptr = D.offset(i as isize);
            *ptr = mCalloc(size, size_of::<f64>() as u64) as *mut f64;
        }
        for i in 0..size {
            for j in 0..size {
                let row = D.offset(i as isize);
                let r2 = *row;
                let ptr = r2.offset(j as isize);
                *ptr = dm[[i as usize, j as usize]];
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
        let mut nniCount : i32 = 0;
        let mut sprCount : i32 = 0;
        t = ImproveTree(&mut options, t, D, A, &mut nniCount, &mut sprCount, options.fpO_stat_file);
        let mut tree_output = vec![0u8; (size << 10) as usize]; //Vec::<u8>::with_capacity((size << 10) as usize);
        NewickPrintTreeStr(t, tree_output.as_mut_ptr() as *mut i8, 2);
        let result = CString::from_vec_unchecked(tree_output);
        return result.to_str().unwrap().to_owned();
    }
}

pub fn parse_newick(taxon_set: &mut TaxonSet, newick: &str) -> Tree {
    let mut taxa : Vec<i32> = vec![-42];
    let mut parents : Vec<i32> = vec![0];
    let mut support : Vec<f64> = vec![-1.0];
    let mut lengths : Vec<f64> = vec![-1.0];
    let mut childcount : Vec<u32> = vec![0];
    let mut firstchild : Vec<i32> = vec![-1];
    let mut nextsib : Vec<i32> = vec![-1];
    // we just reuse TreeSwift's logic
    let mut n : usize = 0; // the current node
    let mut chars = newick.chars().fuse().peekable();
    // let mut parse_length = false;
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
            let mut ls  = "".to_string();
            loop {
                match chars.peek() {
                    Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    },
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
                    },
                    Some(_) => {
                        ts.push(chars.next().unwrap());
                    }
                }
            }
            if childcount[n] == 0 {
                taxa[n] = taxon_set.request(ts) as i32;
                support[n] = 1.0;
            } else {
                support[n] = ts.parse::<f64>().unwrap();
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
    let tree = Tree { taxa, parents, support, lengths, firstchild, nextsib, childcount, fake_root };
    return tree;
}

#[cfg(test)]
mod tests {
    use super::*;
    
 
}

fn main() {
    println!("Hello, world!");
}
