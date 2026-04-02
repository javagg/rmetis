//! K-way Fiduccia-Mattheyses refinement.
//!
//! Generalizes the 2-way FM algorithm to k partitions. For each boundary
//! vertex, computes the gain of moving to each neighboring partition and
//! selects the best feasible move.

use crate::graph::Graph;
use crate::types::Idx;
use priority_queue::PriorityQueue;

/// FM refinement for a k-way partition. Returns final edge cut.
pub fn fm_kway_refine(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],  // [nparts * ncon]
    ubvec: &[f64],   // [ncon]
    niter: usize,
) -> Idx {
    let ncon = graph.ncon;
    let mut pw = graph.partition_weights(part, nparts);

    for _iter in 0..niter {
        if !fm_kway_pass(graph, part, nparts, &mut pw, tpwgts, ubvec, ncon) {
            break;
        }
    }

    graph.edge_cut(part)
}

/// One FM pass for k-way. Returns true if partition improved.
fn fm_kway_pass(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    pw: &mut Vec<Vec<Idx>>,
    tpwgts: &[f64],
    ubvec: &[f64],
    ncon: usize,
) -> bool {
    let n = graph.nvtxs;

    // For each boundary vertex, find best partition to move to and gain
    let mut best_moves: Vec<Option<(usize, Idx, Idx)>> = vec![None; n]; // (target_part, gain, pq_key)

    let compute_best_move = |v: usize, part: &[Idx]| -> Option<(usize, Idx)> {
        let pv = part[v] as usize;
        let mut best_p = usize::MAX;
        let mut best_g = 0i32;

        // Count edges to each neighboring partition
        let s = graph.xadj[v] as usize;
        let e = graph.xadj[v + 1] as usize;
        let mut edge_to = vec![0i32; nparts];
        for j in s..e {
            let u = graph.adjncy[j] as usize;
            let w = graph.edge_weight_at(j);
            edge_to[part[u] as usize] += w;
        }

        for p in 0..nparts {
            if p == pv { continue; }
            if edge_to[p] == 0 { continue; } // no direct connection
            let gain = edge_to[p] - edge_to[pv];
            if gain > best_g {
                best_g = gain;
                best_p = p;
            }
        }
        if best_p == usize::MAX { None } else { Some((best_p, best_g)) }
    };

    // Build priority queue of boundary vertices
    let mut pq: PriorityQueue<usize, Idx> = PriorityQueue::new();
    for v in 0..n {
        let pv = part[v];
        let on_boundary = graph.neighbors(v).iter().any(|&u| part[u as usize] != pv);
        if on_boundary {
            if let Some((tp, gain)) = compute_best_move(v, part) {
                best_moves[v] = Some((v, tp as Idx, gain));
                pq.push(v, gain);
            }
        }
    }

    let mut locked = vec![false; n];
    let mut moves: Vec<(usize, Idx, Idx)> = Vec::new();
    let mut cut = graph.edge_cut(part);
    let mut best_cut = cut;
    let mut best_moves_len = 0usize;

    while let Some((v, gain)) = pq.pop() {
        if locked[v] || gain <= 0 { break; }

        let bm = match best_moves[v] {
            Some(bm) => bm,
            None => continue,
        };
        let to = bm.1 as usize;
        let from = part[v] as usize;
        let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();

        // Balance feasibility
        let feasible = (0..ncon).all(|c| {
            let new_w = pw[to][c] + vw[c];
            new_w as f64 <= tpwgts[to * ncon + c] * ubvec[c]
        });
        if !feasible { continue; }

        // Perform move
        part[v] = to as Idx;
        for c in 0..ncon {
            pw[from][c] -= vw[c];
            pw[to][c] += vw[c];
        }
        locked[v] = true;
        moves.push((v, from as Idx, to as Idx));
        cut -= gain;

        if cut < best_cut {
            best_cut = cut;
            best_moves_len = moves.len();
        }

        // Update neighbors
        let s = graph.xadj[v] as usize;
        let e = graph.xadj[v + 1] as usize;
        for j in s..e {
            let u = graph.adjncy[j] as usize;
            if locked[u] { continue; }
            if let Some((tp, g)) = compute_best_move(u, part) {
                best_moves[u] = Some((u, tp as Idx, g));
                if pq.get(&u).is_some() {
                    pq.change_priority(&u, g);
                } else {
                    pq.push(u, g);
                }
            } else {
                pq.remove(&u);
                best_moves[u] = None;
            }
        }
    }

    // Rollback moves beyond best
    for &(v, from, _to) in moves[best_moves_len..].iter().rev() {
        let cur = part[v] as usize;
        let vw: Vec<Idx> = (0..ncon).map(|c| graph.vwgt_at(v, c)).collect();
        for c in 0..ncon {
            pw[cur][c] -= vw[c];
            pw[from as usize][c] += vw[c];
        }
        part[v] = from;
    }

    best_moves_len > 0 && best_cut < graph.edge_cut(part) + 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    fn grid3x3() -> Graph {
        // 3x3 grid, 9 vertices
        // 0-1-2
        // 3-4-5
        // 6-7-8
        let xadj = vec![0, 2, 4, 6, 9, 12, 14, 16, 19, 21];
        let adjncy = vec![
            1,3,       // 0
            0,2,4,     // 1 (but degree 3, not 2 -- fix)
            1,5,       // 2
            0,4,6,     // 3
            1,3,5,7,   // 4
            2,4,8,     // 5
            3,7,       // 6
            4,6,8,     // 7 (degree 3)
            5,7,       // 8
        ];
        // Rebuild correct adjacency
        let edges = vec![
            (0,1),(0,3),
            (1,2),(1,4),
            (2,5),
            (3,4),(3,6),
            (4,5),(4,7),
            (5,8),
            (6,7),
            (7,8),
        ];
        let n = 9;
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (u, v) in &edges {
            adj[*u].push(*v);
            adj[*v].push(*u);
        }
        let mut xadj2 = vec![0i32; n + 1];
        let mut adjncy2 = Vec::new();
        for v in 0..n {
            adj[v].sort();
            xadj2[v + 1] = xadj2[v] + adj[v].len() as i32;
            for &u in &adj[v] { adjncy2.push(u as i32); }
        }
        Graph::new_unweighted(n, xadj2, adjncy2).unwrap()
    }

    #[test]
    fn kway_does_not_increase_cut() {
        let g = grid3x3();
        let nparts = 3;
        let tpwgts: Vec<f64> = vec![3.0; nparts];
        let ubvec = vec![1.1];
        let mut part: Vec<i32> = (0..9).map(|v| (v / 3) as i32).collect();
        let cut_before = g.edge_cut(&part);
        let cut_after = fm_kway_refine(&g, &mut part, nparts, &tpwgts, &ubvec, 5);
        assert!(cut_after <= cut_before, "{} -> {}", cut_before, cut_after);
    }
}
