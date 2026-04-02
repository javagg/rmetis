//! Fiduccia-Mattheyses (FM) refinement for 2-way partitioning.
//!
//! Classic FM with priority-queue-based vertex selection and rollback to
//! best cut seen during each pass.

use crate::graph::Graph;
use crate::refine::greedy::compute_move_gain;
use crate::types::Idx;
use priority_queue::PriorityQueue;

/// FM refinement for a 2-way partition. Returns final edge cut.
pub fn fm_refine(
    graph: &Graph,
    part: &mut Vec<Idx>,
    tpwgts: &[f64],  // [2 * ncon]
    ubvec: &[f64],   // [ncon]
    niter: usize,
) -> Idx {
    let ncon = graph.ncon;
    let mut pw = graph.partition_weights(part, 2); // [2][ncon]

    for _iter in 0..niter {
        if !fm_pass(graph, part, &mut pw, tpwgts, ubvec, ncon) {
            break;
        }
    }

    graph.edge_cut(part)
}

/// Execute one FM pass. Returns true if the partition improved.
fn fm_pass(
    graph: &Graph,
    part: &mut Vec<Idx>,
    pw: &mut Vec<Vec<Idx>>,
    tpwgts: &[f64],
    ubvec: &[f64],
    ncon: usize,
) -> bool {
    let n = graph.nvtxs;

    // Compute initial gains for all boundary vertices
    let mut gains: Vec<Idx> = vec![0; n];
    for v in 0..n {
        let pv = part[v];
        let other = 1 - pv;
        gains[v] = compute_move_gain(graph, v, part, pv, other);
    }

    // Priority queue keyed by gain (only boundary vertices)
    let mut pq: PriorityQueue<usize, Idx> = PriorityQueue::new();
    for v in 0..n {
        if graph.neighbors(v).iter().any(|&u| part[u as usize] != part[v]) {
            pq.push(v, gains[v]);
        }
    }

    let mut locked = vec![false; n];
    let mut moves: Vec<(usize, Idx, Idx)> = Vec::new(); // (vertex, from, to)
    let mut cut = graph.edge_cut(part);
    let mut best_cut = cut;
    let mut best_moves_len = 0usize;
    let mut cumulative_gain = 0i32;

    while let Some((v, gain)) = pq.pop() {
        if locked[v] {
            continue;
        }

        let from = part[v];
        let to = 1 - from;
        let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();

        // Balance check
        let feasible = (0..ncon).all(|c| {
            let new_w = pw[to as usize][c] + vw[c];
            new_w as f64 <= tpwgts[to as usize * ncon + c] * ubvec[c]
        });
        if !feasible {
            continue;
        }

        // Perform move
        part[v] = to;
        for c in 0..ncon {
            pw[from as usize][c] -= vw[c];
            pw[to as usize][c] += vw[c];
        }
        locked[v] = true;
        moves.push((v, from, to));
        cumulative_gain += gain;
        cut -= gain;

        if cut < best_cut {
            best_cut = cut;
            best_moves_len = moves.len();
        }

        // Update gains for unlocked neighbors
        let s = graph.xadj[v] as usize;
        let e = graph.xadj[v + 1] as usize;
        for j in s..e {
            let u = graph.adjncy[j] as usize;
            if locked[u] {
                continue;
            }
            let pu = part[u];
            let other_u = 1 - pu;
            gains[u] = compute_move_gain(graph, u, part, pu, other_u);
            if pq.get(&u).is_some() {
                pq.change_priority(&u, gains[u]);
            } else {
                // Add to pq if now on boundary
                if graph.neighbors(u).iter().any(|&w| part[w as usize] != pu) {
                    pq.push(u, gains[u]);
                }
            }
        }
    }

    // Rollback moves beyond best position
    let improved = best_moves_len > 0 && best_cut < graph.edge_cut(part);

    // Undo moves after best_moves_len
    for &(v, from, _to) in moves[best_moves_len..].iter().rev() {
        let cur_part = part[v];
        let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();
        pw[cur_part as usize].iter_mut().zip(vw.iter()).for_each(|(w, &vw)| *w -= vw);
        pw[from as usize].iter_mut().zip(vw.iter()).for_each(|(w, &vw)| *w += vw);
        part[v] = from;
    }

    improved
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    fn cycle6() -> Graph {
        // 0-1-2-3-4-5-0
        Graph::new_unweighted(
            6,
            vec![0, 2, 4, 6, 8, 10, 12],
            vec![1,5, 0,2, 1,3, 2,4, 3,5, 4,0],
        ).unwrap()
    }

    #[test]
    fn fm_does_not_increase_cut() {
        let g = cycle6();
        let mut part = vec![0,0,0,1,1,1];
        let tpwgts = vec![3.0, 3.0];
        let ubvec = vec![1.05];
        let cut_before = g.edge_cut(&part);
        let cut_after = fm_refine(&g, &mut part, &tpwgts, &ubvec, 10);
        assert!(cut_after <= cut_before, "cut increased: {} -> {}", cut_before, cut_after);
    }
}
