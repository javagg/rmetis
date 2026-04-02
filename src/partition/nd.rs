//! Nested dissection ordering (NodeND).
//!
//! Produces a fill-reducing vertex ordering for sparse matrix factorization.
//! Uses recursive vertex separator decomposition.

use crate::graph::Graph;
use crate::separator;
use crate::types::{Idx, MetisError, NDResult, Options};
use rand::rngs::SmallRng;
use rand::SeedableRng;

/// Minimum graph size for direct ordering (skip recursion)
const ND_THRESHOLD: usize = 100;

/// Compute a nested dissection ordering of `graph`.
pub fn node_nd(
    graph: &Graph,
    options: &Options,
) -> Result<NDResult, MetisError> {
    let n = graph.nvtxs;
    let seed = if options.seed == 0 {
        let mut rng = SmallRng::from_entropy();
        use rand::Rng;
        rng.gen::<u64>()
    } else {
        options.seed
    };
    let mut rng = SmallRng::seed_from_u64(seed);

    let mut perm = vec![0i32; n];
    let mut iperm = vec![0i32; n];
    let mut counter = 0usize;

    // The ordering assigns "new positions" in bottom-up order:
    // subgraphs first, then their separator last.
    nd_recurse(graph, &(0..n).collect::<Vec<_>>(), options, &mut rng, &mut perm, &mut counter);

    // Build inverse permutation
    for v in 0..n {
        iperm[perm[v] as usize] = v as Idx;
    }

    Ok(NDResult { perm, iperm })
}

/// Recursive nested dissection on a subgraph.
///
/// `orig_ids[i]` = original vertex index of local vertex i.
/// `counter` tracks the next available position in the new ordering.
fn nd_recurse(
    graph: &Graph,
    orig_ids: &[usize],
    options: &Options,
    rng: &mut SmallRng,
    perm: &mut Vec<Idx>,
    counter: &mut usize,
) {
    let n = graph.nvtxs;

    if n <= ND_THRESHOLD {
        // Base case: assign sequential ordering
        for &orig in orig_ids {
            perm[orig] = *counter as Idx;
            *counter += 1;
        }
        return;
    }

    // Find vertex separator
    let (left_local, right_local, sep_local) =
        separator::find_node_separator(graph, options, rng);

    // Map local indices back to original
    let left_orig: Vec<usize> = left_local.iter().map(|&v| orig_ids[v]).collect();
    let right_orig: Vec<usize> = right_local.iter().map(|&v| orig_ids[v]).collect();
    let sep_orig: Vec<usize> = sep_local.iter().map(|&v| orig_ids[v]).collect();

    // If separator is too large (degenerate), fall back to sequential
    if sep_local.len() >= n / 2 {
        for &orig in orig_ids {
            perm[orig] = *counter as Idx;
            *counter += 1;
        }
        return;
    }

    // Extract subgraphs and recurse on left, then right
    let part_for_sub: Vec<Idx> = (0..n).map(|v| {
        if left_local.contains(&v) { 0 }
        else if right_local.contains(&v) { 1 }
        else { 2 } // separator
    }).collect();

    if !left_local.is_empty() {
        let (sub_left, _, n2o) = graph.extract_subgraph(&part_for_sub, 0);
        let sub_orig: Vec<usize> = n2o.iter().map(|&v| orig_ids[v]).collect();
        nd_recurse(&sub_left, &sub_orig, options, rng, perm, counter);
    }

    if !right_local.is_empty() {
        let (sub_right, _, n2o) = graph.extract_subgraph(&part_for_sub, 1);
        let sub_orig: Vec<usize> = n2o.iter().map(|&v| orig_ids[v]).collect();
        nd_recurse(&sub_right, &sub_orig, options, rng, perm, counter);
    }

    // Assign separator vertices last
    for &orig in &sep_orig {
        perm[orig] = *counter as Idx;
        *counter += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    fn path10() -> Graph {
        // Path: 0-1-2-...-9
        let n = 10;
        let mut xadj = vec![0i32; n + 1];
        let mut adjncy = Vec::new();
        for v in 0..n {
            if v > 0 { adjncy.push((v - 1) as i32); }
            if v < n - 1 { adjncy.push((v + 1) as i32); }
            xadj[v + 1] = adjncy.len() as i32;
        }
        Graph::new_unweighted(n, xadj, adjncy).unwrap()
    }

    #[test]
    fn perm_is_valid_permutation() {
        let g = path10();
        let opts = Options { seed: 42, ..Default::default() };
        let result = node_nd(&g, &opts).unwrap();
        let n = g.nvtxs;

        // perm is a permutation of 0..n
        let mut seen = vec![false; n];
        for &p in &result.perm {
            assert!(p >= 0 && (p as usize) < n, "perm value out of range: {}", p);
            assert!(!seen[p as usize], "duplicate perm value: {}", p);
            seen[p as usize] = true;
        }

        // iperm is inverse of perm
        for v in 0..n {
            assert_eq!(result.perm[result.iperm[v] as usize] as usize, v);
        }
    }
}
