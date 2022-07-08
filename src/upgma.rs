use fixedbitset::FixedBitSet;
use ndarray::{Array, Ix2, ShapeBuilder};
use ordered_float::NotNan;
use std::{cmp::Reverse, collections::BinaryHeap};
use tracing::warn;

use crate::tree::Tree;

type MinNotNan = Reverse<NotNan<f64>>;
/// UPGMA*, see Pranjal's thesis section 5.2.1
pub fn upgma_star(distance: &Array<f64, Ix2>, mask: &Array<u32, Ix2>) -> anyhow::Result<Tree> {
    let n = distance.shape()[0];
    let mut m = Array::<f64, _>::zeros((n * 2, n * 2).f());
    let mut known = Array::<u8, _>::zeros((n * 2, n * 2).f());
    let mut pq = BinaryHeap::<(MinNotNan, usize, usize)>::new();
    let mut sizes = vec![1usize; n * 2];
    let mut next_taxa = n; // the taxa to be added next upon join
    let mut absorbed = FixedBitSet::with_capacity(n * 2);
    let mut tree = Tree::rooted_star(n);
    tree.root = 2 * n - 2;
    for i in 0..n - 1 {
        for j in i + 1..n {
            if mask[[i, j]] > 0 {
                m[[i, j]] = distance[[i, j]];
                known[[i, j]] = 1;
                pq.push((Reverse(NotNan::new(m[[i, j]])?), i, j));
            }
        }
    }

    loop {
        let (_d, u, v) = pq.pop().unwrap_or_else(|| {
            warn!("UPGMA*: no more pairs to join, randomly joining pairs");
            let to_join: Vec<usize> = (0..next_taxa)
                .filter(|it| !absorbed.contains(*it))
                .take(2)
                .collect();
            (Reverse(NotNan::new(0.0).unwrap()), to_join[0], to_join[1])
        });
        if absorbed[u] || absorbed[v] {
            continue;
        }
        let new_taxon = next_taxa;
        next_taxa += 1;
        // mark u and v as dead
        absorbed.insert(u);
        absorbed.insert(v);
        sizes[new_taxon] = sizes[u] + sizes[v];
        // maintain the tree
        tree.childcount[new_taxon] = 2;
        tree.firstchild[new_taxon] = u as i32;
        tree.nextsib[u] = v as i32;
        // update distances
        for i in 0..new_taxon {
            if absorbed[i] {
                continue;
            }
            let (u_, i1_) = (std::cmp::min(u, i), std::cmp::max(u, i));
            let (v_, i2_) = (std::cmp::min(v, i), std::cmp::max(v, i));
            match (known[[u_, i1_]], known[[v_, i2_]]) {
                (1, 1) => {
                    let upper = sizes[u] as f64 * m[[u_, i1_]] + sizes[v] as f64 * m[[v_, i2_]];
                    m[[i, new_taxon]] = upper / sizes[new_taxon] as f64
                }
                (1, 0) => m[[i, new_taxon]] = m[[u_, i1_]],
                (0, 1) => m[[i, new_taxon]] = m[[v_, i2_]],
                (0, 0) => {
                    continue;
                }
                (_, _) => {
                    panic!("UPGMA*: corrupted state");
                }
            };
            known[[i, new_taxon]] = 1;
            pq.push((Reverse(NotNan::new(m[[i, new_taxon]])?), i, new_taxon));
        }
        tree.lengths[u] = m[[u, new_taxon]];
        tree.lengths[v] = m[[v, new_taxon]];
        if next_taxa >= 2 * n - 1 {
            let tt = tree.lengths[u] + tree.lengths[v];
            tree.lengths[u] = tt;
            tree.lengths[v] = tt;
            break;
        }
    }
    Ok(tree)
}
