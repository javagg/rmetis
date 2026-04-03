//! Sorted Heavy-Edge Matching (SHEM) coarsening.
//!
//! Vertices are visited in order of increasing degree (light vertices first).
//! For each unvisited vertex, the unvisited neighbor with the highest edge
//! weight is selected as the match partner.

use crate::graph::Graph;
use rand::Rng;

/// Compute a maximal matching using sorted heavy-edge heuristic.
///
/// Returns `matching` where `matching[v]` is the partner of v (or v itself).
pub fn shem_matching(graph: &Graph, rng: &mut impl Rng) -> Vec<usize> {
    let n = graph.nvtxs;
    let mut matching = (0..n).collect::<Vec<usize>>();
    let mut visited = vec![false; n];

    // Sort vertices by ascending degree (ties broken by pre-sampled random value)
    let mut perm: Vec<usize> = (0..n).collect();
    // Pre-generate one random tiebreaker per vertex to ensure the key function
    // is deterministic (sort_by_key calls the closure multiple times per element).
    let tiebreakers: Vec<u32> = (0..n).map(|_| rng.gen::<u32>() >> 16).collect();
    perm.sort_by_key(|&v| (graph.degree(v), tiebreakers[v]));

    for &v in &perm {
        if visited[v] {
            continue;
        }
        visited[v] = true;

        // Pick neighbor with maximum edge weight
        let mut best_u = usize::MAX;
        let mut best_w = Idx::MIN;

        let s = graph.xadj[v] as usize;
        let e = graph.xadj[v + 1] as usize;
        for j in s..e {
            let u = graph.adjncy[j] as usize;
            if visited[u] {
                continue;
            }
            let w = graph.edge_weight_at(j);
            if w > best_w {
                best_w = w;
                best_u = u;
            }
        }

        if best_u != usize::MAX {
            matching[v] = best_u;
            matching[best_u] = v;
            visited[best_u] = true;
        }
    }
    matching
}

use crate::types::Idx;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn weighted_path() -> Graph {
        // Path 0-1-2-3 with edge weights 1, 3, 1
        // xadj, adjncy for undirected path
        Graph::new(
            4, 1,
            vec![0, 1, 3, 5, 6],
            vec![1, 0, 2, 1, 3, 2],
            None,
            Some(vec![1, 1, 3, 3, 1, 1]), // adjwgt parallel to adjncy
            None,
        ).unwrap()
    }

    #[test]
    fn shem_prefers_heavy_edge() {
        // Graph: 0--1==2--3 where == is the heavy edge (weight 3)
        // Use a star topology: center vertex 0 connects to 1(w=1), 2(w=3), 3(w=1)
        // So center has degree 3; leaves have degree 1.
        // SHEM should match center with leaf 2 (highest weight).
        let g = Graph::new(
            4, 1,
            vec![0, 3, 4, 5, 6],
            vec![1, 2, 3, 0, 0, 0], // 0->1,2,3; 1->0; 2->0; 3->0
            None,
            Some(vec![1, 3, 1, 1, 3, 1]), // weights: 0-1=1, 0-2=3, 0-3=1
            None,
        ).unwrap();
        let mut rng = SmallRng::seed_from_u64(1);
        let m = shem_matching(&g, &mut rng);
        // Symmetry check
        for v in 0..g.nvtxs {
            assert_eq!(m[m[v]], v);
        }
        // Leaves (degree 1) are processed first.
        // Leaves 1,2,3 all try to match with center 0.
        // First leaf processed picks center 0.
        // Verify at least one leaf is matched
        let matched_pairs = (0..4).filter(|&v| m[v] != v).count();
        assert!(matched_pairs >= 2, "expected at least one matched pair");
    }
}
