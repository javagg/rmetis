//! Greedy boundary refinement.
//!
//! Iteratively moves boundary vertices to the neighboring partition that
//! gives the highest gain, as long as the move satisfies the balance
//! constraint. Repeats until no improvement or `niter` passes are done.

use crate::comm::{Comm, SingleComm};
use crate::graph::Graph;
use crate::types::Idx;

/// Greedy k-way refinement. Returns final edge cut.
pub fn greedy_refine(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],  // [nparts * ncon]
    ubvec: &[f64],   // [ncon]
    niter: usize,
) -> Idx {
    let single = SingleComm;
    greedy_refine_with_comm(graph, part, nparts, tpwgts, ubvec, niter, &single)
}

/// Greedy k-way refinement with MPI parallelism.
///
/// Boundary vertices are distributed across ranks (each rank owns vertices
/// `v` where `v % nprocs == rank`). Each rank proposes moves for its subset;
/// after each pass all ranks exchange proposed moves via `all_gather` and
/// apply them atomically, then re-sync partition weights.
///
/// A balance-repair pass runs first: if any partition exceeds its weight
/// budget, vertices are moved to the lightest feasible neighbor even if
/// the move does not improve the cut.
///
/// Returns the global edge cut (identical on all ranks).
pub fn greedy_refine_with_comm(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],
    ubvec: &[f64],
    niter: usize,
    comm: &dyn Comm,
) -> Idx {
    let ncon = graph.ncon;
    let nprocs = comm.size().max(1) as usize;
    let rank = comm.rank() as usize;

    // Current partition weights — all ranks start with the same view
    let mut pw = graph.partition_weights(part, nparts); // [nparts][ncon]

    // Flat pw for all_reduce sync: layout [p0c0, p0c1, ..., p1c0, ...]
    let pw_len = nparts * ncon;

    // ── Balance-repair pass ──────────────────────────────────────────────────
    // Move boundary vertices away from overweight partitions to feasible
    // neighbours even if the move has zero or negative cut gain. This runs
    // before the cut-optimisation loop so that the main loop starts with a
    // balanced partition.
    {
        let boundary = graph.boundary_vertices(part);
        let mut local_moves: Vec<Idx> = Vec::new();

        for v in boundary.iter().copied().filter(|&v| v % nprocs == rank) {
            let pv = part[v] as usize;
            let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();

            // Only attempt to move if the source partition is overweight
            let source_overweight = (0..ncon).any(|c| {
                pw[pv][c] as f64 > tpwgts[pv * ncon + c] * ubvec[c]
            });
            if !source_overweight {
                continue;
            }

            // Among feasible destination partitions, pick the one with best gain
            // (or least-negative gain if no positive option exists)
            let mut best_gain = Idx::MIN;
            let mut best_p = usize::MAX;

            for p in 0..nparts {
                if p == pv {
                    continue;
                }
                let feasible = (0..ncon).all(|c| {
                    let new_w = pw[p][c] + vw[c];
                    new_w as f64 <= tpwgts[p * ncon + c] * ubvec[c]
                });
                if !feasible {
                    continue;
                }
                let gain = compute_move_gain(graph, v, part, pv as Idx, p as Idx);
                if gain > best_gain {
                    best_gain = gain;
                    best_p = p;
                }
            }

            if best_p != usize::MAX {
                local_moves.push(v as Idx);
                local_moves.push(best_p as Idx);
                for c in 0..ncon {
                    pw[pv][c] -= vw[c];
                    pw[best_p][c] += vw[c];
                }
            }
        }

        // Exchange balance-repair moves
        let all_moves = comm.all_gather_i32_vec(&local_moves);
        for rank_moves in all_moves {
            let mut i = 0;
            while i + 1 < rank_moves.len() {
                let v = rank_moves[i] as usize;
                let new_p = rank_moves[i + 1] as usize;
                part[v] = new_p as Idx;
                i += 2;
            }
        }
        if nprocs > 1 {
            pw = graph.partition_weights(part, nparts);
        }
    }

    // ── Cut-optimisation passes ──────────────────────────────────────────────
    for _iter in 0..niter {
        let boundary = graph.boundary_vertices(part);

        let mut local_moves: Vec<Idx> = Vec::new();

        for v in boundary.iter().copied().filter(|&v| v % nprocs == rank) {
            let pv = part[v] as usize;
            let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();

            let mut best_gain = 0i32;
            let mut best_p = pv;

            for p in 0..nparts {
                if p == pv {
                    continue;
                }
                let feasible = (0..ncon).all(|c| {
                    let new_w = pw[p][c] + vw[c];
                    new_w as f64 <= tpwgts[p * ncon + c] * ubvec[c]
                });
                if !feasible {
                    continue;
                }
                let gain = compute_move_gain(graph, v, part, pv as Idx, p as Idx);
                if gain > best_gain {
                    best_gain = gain;
                    best_p = p;
                }
            }

            if best_p != pv {
                local_moves.push(v as Idx);
                local_moves.push(best_p as Idx);
                for c in 0..ncon {
                    pw[pv][c] -= vw[c];
                    pw[best_p][c] += vw[c];
                }
            }
        }

        let all_moves = comm.all_gather_i32_vec(&local_moves);

        let mut any_move = false;
        for rank_moves in all_moves {
            let mut i = 0;
            while i + 1 < rank_moves.len() {
                let v = rank_moves[i] as usize;
                let new_p = rank_moves[i + 1] as usize;
                let old_p = part[v] as usize;
                if old_p != new_p {
                    part[v] = new_p as Idx;
                    any_move = true;
                }
                i += 2;
            }
        }

        if !any_move {
            break;
        }

        if nprocs > 1 {
            pw = graph.partition_weights(part, nparts);
        }
    }

    graph.edge_cut(part)
}

/// Compute the gain of moving vertex v from `from_part` to `to_part`.
/// gain > 0 means the move reduces the edge cut.
pub fn compute_move_gain(
    graph: &Graph,
    v: usize,
    part: &[Idx],
    from_part: Idx,
    to_part: Idx,
) -> Idx {
    let s = graph.xadj[v] as usize;
    let e = graph.xadj[v + 1] as usize;
    let mut gain = 0i32;
    for j in s..e {
        let u = graph.adjncy[j] as usize;
        let w = graph.edge_weight_at(j);
        let pu = part[u];
        if pu == from_part {
            gain -= w;
        } else if pu == to_part {
            gain += w;
        }
    }
    gain
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    fn cycle4() -> Graph {
        Graph::new_unweighted(4, vec![0, 2, 4, 6, 8], vec![1, 3, 0, 2, 1, 3, 0, 2]).unwrap()
    }

    #[test]
    fn refine_does_not_increase_cut() {
        let g = cycle4();
        let mut part = vec![0, 1, 0, 1];
        let tpwgts = vec![2.0, 2.0];
        let ubvec = vec![1.03];
        let cut_before = g.edge_cut(&part);
        let cut_after = greedy_refine(&g, &mut part, 2, &tpwgts, &ubvec, 10);
        assert!(cut_after <= cut_before);
    }

    #[test]
    fn refine_single_comm_matches_serial() {
        let g = cycle4();
        let mut part_a = vec![0, 1, 0, 1];
        let mut part_b = vec![0, 1, 0, 1];
        let tpwgts = vec![2.0, 2.0];
        let ubvec = vec![1.03];
        let cut_a = greedy_refine(&g, &mut part_a, 2, &tpwgts, &ubvec, 10);
        let cut_b = greedy_refine_with_comm(&g, &mut part_b, 2, &tpwgts, &ubvec, 10, &SingleComm);
        assert_eq!(cut_a, cut_b);
        assert_eq!(part_a, part_b);
    }
}
