//! Graph-growing (BFS-based) initial partitioning.

use crate::graph::Graph;
use crate::types::Idx;
use rand::Rng;
use std::collections::VecDeque;

/// Partition graph into `nparts` using BFS growth from random seeds.
///
/// `tpwgts[p * ncon + c]` = target weight (as float, already scaled to integer
/// total weight) for partition p, constraint c.
pub fn grow_partition(
    graph: &Graph,
    nparts: usize,
    tpwgts: &[f64],    // shape [nparts * ncon], target integer weights
    rng: &mut impl Rng,
) -> Vec<Idx> {
    let n = graph.nvtxs;
    let ncon = graph.ncon;
    let mut part = vec![-1i32; n];
    let mut part_weight = vec![vec![0i64; ncon]; nparts];

    // Pick random seed vertices (distinct)
    let mut seeds: Vec<usize> = Vec::with_capacity(nparts);
    let mut used = vec![false; n];
    while seeds.len() < nparts {
        let v = rng.gen_range(0..n);
        if !used[v] {
            used[v] = true;
            seeds.push(v);
        }
    }

    // BFS queues per partition
    let mut queues: Vec<VecDeque<usize>> = (0..nparts).map(|p| {
        let mut q = VecDeque::new();
        q.push_back(seeds[p]);
        q
    }).collect();

    // Assign seeds
    for (p, &s) in seeds.iter().enumerate() {
        part[s] = p as Idx;
        for c in 0..ncon {
            part_weight[p][c] += graph.vwgt_at(s, c) as i64;
        }
    }

    // BFS growth: round-robin across partitions
    let mut active = nparts;
    while active > 0 {
        active = 0;
        for p in 0..nparts {
            // Check if partition p still needs more weight
            let at_target = (0..ncon).all(|c| {
                part_weight[p][c] as f64 >= tpwgts[p * ncon + c]
            });
            if at_target {
                continue;
            }

            // Try to grow partition p
            let mut grew = false;
            while let Some(v) = queues[p].pop_front() {
                let s = graph.xadj[v] as usize;
                let e = graph.xadj[v + 1] as usize;
                for j in s..e {
                    let u = graph.adjncy[j] as usize;
                    if part[u] == -1 {
                        part[u] = p as Idx;
                        for c in 0..ncon {
                            part_weight[p][c] += graph.vwgt_at(u, c) as i64;
                        }
                        queues[p].push_back(u);
                        grew = true;
                    }
                }
                if grew {
                    break; // process one vertex per round
                }
            }
            if grew {
                active += 1;
            }
        }
    }

    // Assign remaining unassigned vertices (isolated or disconnected)
    for v in 0..n {
        if part[v] == -1 {
            // Find partition with minimum current weight (constraint 0)
            let best_p = (0..nparts)
                .min_by_key(|&p| part_weight[p][0])
                .unwrap_or(0);
            part[v] = best_p as Idx;
            for c in 0..ncon {
                part_weight[best_p][c] += graph.vwgt_at(v, c) as i64;
            }
        }
    }

    part
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn grid4() -> Graph {
        // 2x2 grid: 0-1, 0-2, 1-3, 2-3
        Graph::new_unweighted(
            4,
            vec![0, 2, 4, 6, 8],
            vec![1, 2, 0, 3, 0, 3, 1, 2],
        ).unwrap()
    }

    #[test]
    fn all_vertices_assigned() {
        let g = grid4();
        let tpwgts = vec![2.0, 2.0]; // 2 parts, equal weight
        let mut rng = SmallRng::seed_from_u64(42);
        let part = grow_partition(&g, 2, &tpwgts, &mut rng);
        assert!(part.iter().all(|&p| p >= 0 && p < 2));
    }

    #[test]
    fn partition_count() {
        let g = grid4();
        let tpwgts = vec![2.0, 2.0];
        let mut rng = SmallRng::seed_from_u64(7);
        let part = grow_partition(&g, 2, &tpwgts, &mut rng);
        let c0 = part.iter().filter(|&&p| p == 0).count();
        let c1 = part.iter().filter(|&&p| p == 1).count();
        assert_eq!(c0 + c1, 4);
    }
}
