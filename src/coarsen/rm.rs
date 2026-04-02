//! Random Matching (RM) coarsening.
//!
//! Vertices are visited in random order; the first unvisited neighbor of each
//! unvisited vertex is selected as its match partner.

use crate::graph::Graph;
use rand::Rng;

/// Compute a maximal matching using random ordering.
///
/// Returns a `matching` array where `matching[v]` is the partner of vertex v.
/// Unmatched vertices are self-matched: `matching[v] == v`.
pub fn random_matching(graph: &Graph, rng: &mut impl Rng) -> Vec<usize> {
    let n = graph.nvtxs;
    let mut matching = (0..n).collect::<Vec<usize>>(); // default: self-match
    let mut visited = vec![false; n];

    // Random permutation of vertices
    let mut perm: Vec<usize> = (0..n).collect();
    // Fisher-Yates shuffle
    for i in (1..n).rev() {
        let j = rng.gen_range(0..=i);
        perm.swap(i, j);
    }

    for &v in &perm {
        if visited[v] {
            continue;
        }
        visited[v] = true;
        // Find the first unvisited neighbor
        for &u in graph.neighbors(v) {
            let u = u as usize;
            if !visited[u] {
                matching[v] = u;
                matching[u] = v;
                visited[u] = true;
                break;
            }
        }
    }
    matching
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn path4() -> Graph {
        Graph::new_unweighted(4, vec![0, 1, 3, 5, 6], vec![1, 0, 2, 1, 3, 2]).unwrap()
    }

    #[test]
    fn matching_is_symmetric() {
        let g = path4();
        let mut rng = SmallRng::seed_from_u64(42);
        let m = random_matching(&g, &mut rng);
        for v in 0..g.nvtxs {
            assert_eq!(m[m[v]], v, "matching not symmetric at {}", v);
        }
    }

    #[test]
    fn matching_covers_neighbors() {
        let g = path4();
        let mut rng = SmallRng::seed_from_u64(0);
        let m = random_matching(&g, &mut rng);
        // At least two vertices should be matched (path 0-1-2-3 has 2 edges)
        let matched = m.iter().enumerate().filter(|&(v, &u)| u != v).count();
        assert!(matched >= 2, "expected at least 2 matched vertices, got {}", matched);
    }
}
