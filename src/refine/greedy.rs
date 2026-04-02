//! Greedy boundary refinement.
//!
//! Iteratively moves boundary vertices to the neighboring partition that
//! gives the highest gain, as long as the move satisfies the balance
//! constraint. Repeats until no improvement or `niter` passes are done.

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
    let n = graph.nvtxs;
    let ncon = graph.ncon;

    // Current partition weights
    let mut pw = graph.partition_weights(part, nparts); // [nparts][ncon]

    for _iter in 0..niter {
        let mut improved = false;

        // Collect boundary vertices
        let boundary = graph.boundary_vertices(part);

        for v in boundary {
            let pv = part[v] as usize;
            let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();

            // Compute gain for moving v to each other partition
            let mut best_gain = 0i32;
            let mut best_p = pv;

            for p in 0..nparts {
                if p == pv {
                    continue;
                }

                // Balance check: would partition p be overweight?
                let feasible = (0..ncon).all(|c| {
                    let new_w = pw[p][c] + vw[c];
                    new_w as f64 <= tpwgts[p * ncon + c] * ubvec[c]
                });
                if !feasible {
                    continue;
                }

                // Gain = (edges to p) - (edges to pv)
                let gain = compute_move_gain(graph, v, part, pv as Idx, p as Idx);
                if gain > best_gain {
                    best_gain = gain;
                    best_p = p;
                }
            }

            if best_p != pv {
                // Perform move
                for c in 0..ncon {
                    pw[pv][c] -= vw[c];
                    pw[best_p][c] += vw[c];
                }
                part[v] = best_p as Idx;
                improved = true;
            }
        }

        if !improved {
            break;
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
            gain -= w; // was internal, becomes external
        } else if pu == to_part {
            gain += w; // was external (to to_part), becomes internal
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
        let mut part = vec![0, 1, 0, 1]; // cross-diagonal partition
        let tpwgts = vec![2.0, 2.0];
        let ubvec = vec![1.03];
        let cut_before = g.edge_cut(&part);
        let cut_after = greedy_refine(&g, &mut part, 2, &tpwgts, &ubvec, 10);
        assert!(cut_after <= cut_before);
    }
}
