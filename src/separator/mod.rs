//! Vertex separator computation for nested dissection.

use crate::graph::Graph;
use crate::types::{Idx, Options};
use rand::rngs::SmallRng;

/// Find a vertex separator for the graph using edge-separator-to-vertex-separator conversion.
///
/// Returns `(left, right, separator)` — vertex sets after removing the separator.
pub fn find_node_separator(
    graph: &Graph,
    options: &Options,
    rng: &mut SmallRng,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    use crate::partition::recursive::partition_recursive;
    use crate::types::Real;

    // 1. Bisect the graph
    let tpwgts = [0.5 as Real, 0.5 as Real];
    let ubvec = [1.05f32];
    let result = partition_recursive(graph, 2, Some(&tpwgts), Some(&ubvec), options)
        .unwrap_or_else(|_| {
            // Fallback: trivial bisection
            let n = graph.nvtxs;
            let part: Vec<Idx> = (0..n).map(|v| if v < n / 2 { 0 } else { 1 }).collect();
            crate::types::PartitionResult { part, objval: 0 }
        });

    let part = result.part;

    // 2. Find edge-separator edges: edges crossing the partition boundary
    // Convert to vertex separator: add the endpoint from the smaller side
    let pw = graph.partition_weights(&part, 2);
    let smaller_side = if pw[0][0] <= pw[1][0] { 0i32 } else { 1i32 };

    let mut in_sep = vec![false; graph.nvtxs];

    // For each cut edge, add the endpoint on the smaller side to separator
    for v in 0..graph.nvtxs {
        if part[v] != smaller_side {
            continue;
        }
        let has_cut_neighbor = graph.neighbors(v)
            .iter()
            .any(|&u| part[u as usize] != part[v]);
        if has_cut_neighbor {
            in_sep[v] = true;
        }
    }

    // 3. Build the three sets
    let separator: Vec<usize> = (0..graph.nvtxs).filter(|&v| in_sep[v]).collect();
    let left: Vec<usize> = (0..graph.nvtxs)
        .filter(|&v| !in_sep[v] && part[v] == smaller_side)
        .collect();
    let right: Vec<usize> = (0..graph.nvtxs)
        .filter(|&v| !in_sep[v] && part[v] != smaller_side)
        .collect();

    (left, right, separator)
}
