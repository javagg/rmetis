//! Tests for parallel (MPI comm abstraction) components.
//!
//! These tests exercise the parallel algorithms using `SingleComm` (1 rank),
//! which is the only mode testable without a real MPI runtime. They verify:
//! 1. Comm trait implementations work correctly on SingleComm
//! 2. Parallel algorithms produce the same results as serial when comm is single-rank
//! 3. The full partition_kway pipeline with an explicit comm produces valid results

use rmetis::comm::{Comm, SingleComm};
use rmetis::coarsen;
use rmetis::graph::Graph;
use rmetis::refine::greedy::{greedy_refine, greedy_refine_with_comm};
use rmetis::partition::kway::partition_kway;
use rmetis::{Options, PartitionResult};

fn cycle4() -> Graph {
    Graph::new_unweighted(4, vec![0, 2, 4, 6, 8], vec![1, 3, 0, 2, 1, 3, 0, 2]).unwrap()
}

fn grid_graph(rows: usize, cols: usize) -> Graph {
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjncy = Vec::new();

    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            let mut neighbors = Vec::new();
            if r > 0 { neighbors.push((r - 1) * cols + c); }
            if r + 1 < rows { neighbors.push((r + 1) * cols + c); }
            if c > 0 { neighbors.push(r * cols + c - 1); }
            if c + 1 < cols { neighbors.push(r * cols + c + 1); }
            xadj[v + 1] = xadj[v] + neighbors.len() as i32;
            adjncy.extend(neighbors.iter().map(|&x| x as i32));
        }
    }
    Graph::new_unweighted(n, xadj, adjncy).unwrap()
}

// ─── SingleComm trait tests ──────────────────────────────────────────────────

#[test]
fn single_comm_rank_and_size() {
    let c = SingleComm;
    assert_eq!(c.rank(), 0);
    assert_eq!(c.size(), 1);
}

#[test]
fn single_comm_all_reduce_min_is_identity() {
    let c = SingleComm;
    assert_eq!(c.all_reduce_min_i32(42), 42);
    assert_eq!(c.all_reduce_min_i32(-5), -5);
}

#[test]
fn single_comm_all_reduce_sum_slice_is_copy() {
    let c = SingleComm;
    let local = vec![1, 2, 3, 4];
    let mut out = vec![0; 4];
    c.all_reduce_sum_i32_slice(&local, &mut out);
    assert_eq!(out, local);
}

#[test]
fn single_comm_gather_returns_single_row() {
    let c = SingleComm;
    let result = c.gather_i32_vec(0, &[10, 20, 30]);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], vec![10, 20, 30]);
}

#[test]
fn single_comm_all_gather_returns_single_row() {
    let c = SingleComm;
    let result = c.all_gather_i32_vec(&[5, 6]);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], vec![5, 6]);
}

#[test]
fn single_comm_broadcast_is_noop() {
    let c = SingleComm;
    let mut data = vec![100i32, 200];
    c.broadcast_i32_vec(0, &mut data);
    assert_eq!(data, vec![100, 200]);
}

// ─── Distributed greedy refinement tests ────────────────────────────────────

#[test]
fn greedy_refine_with_single_comm_matches_serial() {
    let g = cycle4();

    let test_cases = [
        vec![0i32, 0, 1, 1],
        vec![0i32, 1, 0, 1],
        vec![0i32, 1, 1, 0],
    ];

    for initial_part in &test_cases {
        let tpwgts = vec![2.0f64, 2.0];
        let ubvec = vec![1.03f64];

        let mut part_serial = initial_part.clone();
        let mut part_parallel = initial_part.clone();

        let cut_serial = greedy_refine(&g, &mut part_serial, 2, &tpwgts, &ubvec, 10);
        let cut_parallel = greedy_refine_with_comm(&g, &mut part_parallel, 2, &tpwgts, &ubvec, 10, &SingleComm);

        assert_eq!(cut_serial, cut_parallel,
            "serial vs parallel cut mismatch for initial {:?}", initial_part);
        assert_eq!(part_serial, part_parallel,
            "serial vs parallel part mismatch for initial {:?}", initial_part);
    }
}

#[test]
fn greedy_refine_with_comm_does_not_increase_cut() {
    let g = grid_graph(4, 4);
    let nparts = 4;
    let tpwgts: Vec<f64> = vec![4.0; nparts]; // uniform target
    let ubvec = vec![1.05f64];

    // A fairly bad initial partition: split into 4 rows
    let mut part: Vec<i32> = (0..16).map(|v| (v / 4) as i32).collect();
    let cut_before = g.edge_cut(&part);
    let cut_after = greedy_refine_with_comm(&g, &mut part, nparts, &tpwgts, &ubvec, 10, &SingleComm);

    assert!(cut_after <= cut_before,
        "greedy_refine_with_comm increased cut: {} -> {}", cut_before, cut_after);
}

// ─── Distributed coarsening tests ───────────────────────────────────────────

#[test]
fn coarsen_with_single_comm_produces_valid_hierarchy() {
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    let g = grid_graph(4, 4);
    let opts = Options { seed: 42, ..Options::default() };
    let mut rng = SmallRng::seed_from_u64(42);

    let levels = coarsen::coarsen_with_comm(&g, 2, &opts, &mut rng, &SingleComm);
    assert!(!levels.is_empty(), "coarsening produced no levels");

    // Each level's cmap must be valid (values in [0, coarse_nvtxs))
    let mut prev_n = g.nvtxs;
    for (i, level) in levels.iter().enumerate() {
        assert_eq!(level.cmap.len(), prev_n,
            "level {} cmap length {} != prev_n {}", i, level.cmap.len(), prev_n);
        let coarse_n = level.graph.nvtxs;
        for &cv in &level.cmap {
            assert!(cv < coarse_n, "level {} cmap value {} out of range (coarse_n={})", i, cv, coarse_n);
        }
        prev_n = coarse_n;
    }
}

// ─── Full pipeline tests ─────────────────────────────────────────────────────

#[test]
fn partition_kway_with_explicit_single_comm_is_valid() {
    let g = grid_graph(4, 4);
    let opts = Options { seed: 1, ncuts: 3, ..Options::for_kway() };

    let result = partition_kway(&g, 4, None, None, &opts, Some(&SingleComm)).unwrap();

    // All vertices assigned to valid partitions
    assert_eq!(result.part.len(), g.nvtxs);
    for &p in &result.part {
        assert!(p >= 0 && p < 4, "invalid partition id {}", p);
    }

    // objval matches computed edge cut
    assert_eq!(result.objval, g.edge_cut(&result.part));

    // Each partition is non-empty
    let mut counts = vec![0usize; 4];
    for &p in &result.part {
        counts[p as usize] += 1;
    }
    for (i, &c) in counts.iter().enumerate() {
        assert!(c > 0, "partition {} is empty", i);
    }
}

#[test]
fn partition_kway_with_none_comm_matches_explicit_single_comm() {
    let g = grid_graph(4, 4);
    let opts = Options { seed: 7, ncuts: 2, ..Options::for_kway() };

    let result_none = partition_kway(&g, 4, None, None, &opts, None).unwrap();
    let result_comm = partition_kway(&g, 4, None, None, &opts, Some(&SingleComm)).unwrap();

    // With the same seed and ncuts=2, both should produce the same result
    assert_eq!(result_none.objval, result_comm.objval,
        "objval differs: none={} comm={}", result_none.objval, result_comm.objval);
}

#[test]
fn partition_kway_parallel_ncuts_covers_all_attempts() {
    // With ncuts=4 and a single-rank comm, rank 0 handles attempts 0,1,2,3.
    // With ncuts=1, only attempt 0 runs. The ncuts=4 result should be <= ncuts=1.
    let g = grid_graph(6, 6);
    let opts_1cut = Options { seed: 42, ncuts: 1, ..Options::for_kway() };
    let opts_4cut = Options { seed: 42, ncuts: 4, ..Options::for_kway() };

    let r1 = partition_kway(&g, 4, None, None, &opts_1cut, Some(&SingleComm)).unwrap();
    let r4 = partition_kway(&g, 4, None, None, &opts_4cut, Some(&SingleComm)).unwrap();

    assert!(r4.objval <= r1.objval,
        "more cuts should not worsen result: 1cut={} 4cut={}", r1.objval, r4.objval);
}
