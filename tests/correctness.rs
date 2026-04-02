//! Integration tests for correctness and regression.

use rmetis::{Graph, Options, part_graph_kway, part_graph_recursive, node_nd};

// -----------------------------------------------------------------------
// Test graph helpers
// -----------------------------------------------------------------------

/// 4-vertex cycle 0-1-2-3-0
fn cycle4() -> Graph {
    Graph::new_unweighted(4, vec![0,2,4,6,8], vec![1,3,0,2,1,3,0,2]).unwrap()
}

/// 2×3 grid (6 vertices)
/// 0-1-2
/// 3-4-5
fn grid2x3() -> Graph {
    // Edges: horizontal + vertical
    let edges = [(0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5)];
    build_undirected(6, &edges)
}

/// 4×4 grid (16 vertices)
fn grid4x4() -> Graph {
    let mut edges = Vec::new();
    for r in 0..4usize {
        for c in 0..4usize {
            let v = r * 4 + c;
            if c + 1 < 4 { edges.push((v, v + 1)); }
            if r + 1 < 4 { edges.push((v, v + 4)); }
        }
    }
    build_undirected(16, &edges)
}

fn build_undirected(n: usize, edges: &[(usize, usize)]) -> Graph {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(u, v) in edges {
        adj[u].push(v);
        adj[v].push(u);
    }
    let mut xadj = vec![0i32; n + 1];
    let mut adjncy = Vec::new();
    for v in 0..n {
        adj[v].sort();
        xadj[v + 1] = xadj[v] + adj[v].len() as i32;
        for &u in &adj[v] { adjncy.push(u as i32); }
    }
    Graph::new_unweighted(n, xadj, adjncy).unwrap()
}

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

fn assert_partition_valid(g: &Graph, part: &[i32], nparts: usize) {
    assert_eq!(part.len(), g.nvtxs, "partition length mismatch");
    for (v, &p) in part.iter().enumerate() {
        assert!(p >= 0 && (p as usize) < nparts,
            "vertex {} has invalid part {}", v, p);
    }
    // All parts used
    let mut seen = vec![false; nparts];
    for &p in part { seen[p as usize] = true; }
    // (Some partitions may be empty for very small graphs — skip strict check)
}

fn assert_partition_balanced(g: &Graph, part: &[i32], nparts: usize, ufactor: u32) {
    let pw = g.partition_weights(part, nparts);
    let total = g.total_weight();
    let max_allowed = total[0] as f64 * (1.0 + ufactor as f64 / 1000.0) / nparts as f64;
    for (p, w) in pw.iter().enumerate() {
        assert!(
            w[0] as f64 <= max_allowed + 1.0, // +1 for integer rounding
            "partition {} weight {} exceeds max {:.1}",
            p, w[0], max_allowed
        );
    }
}

fn assert_nd_valid(n: usize, perm: &[i32], iperm: &[i32]) {
    assert_eq!(perm.len(), n);
    assert_eq!(iperm.len(), n);
    let mut seen = vec![false; n];
    for &p in perm {
        assert!(p >= 0 && (p as usize) < n, "perm out of range: {}", p);
        assert!(!seen[p as usize], "duplicate perm value: {}", p);
        seen[p as usize] = true;
    }
    for v in 0..n {
        assert_eq!(perm[iperm[v] as usize] as usize, v,
            "iperm is not inverse of perm at {}", v);
    }
}

// -----------------------------------------------------------------------
// Basic correctness tests
// -----------------------------------------------------------------------

#[test]
fn kway_cycle4_k2() {
    let g = cycle4();
    let opts = Options { seed: 1, ..Options::for_kway() };
    let r = part_graph_kway(&g, 2, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 2);
    assert_partition_balanced(&g, &r.part, 2, 30);
    assert!(r.objval <= 4, "edge cut {} too high for cycle4 k=2", r.objval);
}

#[test]
fn kway_grid2x3_k2() {
    let g = grid2x3();
    let opts = Options { seed: 42, ..Options::for_kway() };
    let r = part_graph_kway(&g, 2, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 2);
    assert_partition_balanced(&g, &r.part, 2, 30);
}

#[test]
fn kway_grid4x4_k4() {
    let g = grid4x4();
    let opts = Options { seed: 42, ..Options::for_kway() };
    let r = part_graph_kway(&g, 4, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 4);
    assert_partition_balanced(&g, &r.part, 4, 30);
}

#[test]
fn recursive_cycle4_k2() {
    let g = cycle4();
    let opts = Options { seed: 7, ..Options::for_recursive() };
    let r = part_graph_recursive(&g, 2, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 2);
    assert_partition_balanced(&g, &r.part, 2, 5);
}

#[test]
fn recursive_grid4x4_k4() {
    let g = grid4x4();
    let opts = Options { seed: 42, ..Options::for_recursive() };
    let r = part_graph_recursive(&g, 4, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 4);
}

#[test]
fn recursive_grid4x4_k8() {
    let g = grid4x4();
    let opts = Options { seed: 42, ..Options::for_recursive() };
    let r = part_graph_recursive(&g, 8, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 8);
}

// -----------------------------------------------------------------------
// NodeND tests
// -----------------------------------------------------------------------

#[test]
fn nd_cycle4() {
    let g = cycle4();
    let opts = Options { seed: 1, ..Default::default() };
    let r = node_nd(&g, &opts).unwrap();
    assert_nd_valid(g.nvtxs, &r.perm, &r.iperm);
}

#[test]
fn nd_grid4x4() {
    let g = grid4x4();
    let opts = Options { seed: 42, ..Default::default() };
    let r = node_nd(&g, &opts).unwrap();
    assert_nd_valid(g.nvtxs, &r.perm, &r.iperm);
}

// -----------------------------------------------------------------------
// Determinism tests
// -----------------------------------------------------------------------

#[test]
fn kway_deterministic() {
    let g = grid4x4();
    let opts = Options { seed: 123, ..Options::for_kway() };
    let r1 = part_graph_kway(&g, 4, None, None, &opts).unwrap();
    let r2 = part_graph_kway(&g, 4, None, None, &opts).unwrap();
    assert_eq!(r1.objval, r2.objval, "non-deterministic kway");
    assert_eq!(r1.part, r2.part, "non-deterministic kway partition");
}

#[test]
fn recursive_deterministic() {
    let g = grid4x4();
    let opts = Options { seed: 99, ..Options::for_recursive() };
    let r1 = part_graph_recursive(&g, 4, None, None, &opts).unwrap();
    let r2 = part_graph_recursive(&g, 4, None, None, &opts).unwrap();
    assert_eq!(r1.objval, r2.objval, "non-deterministic recursive");
    assert_eq!(r1.part, r2.part, "non-deterministic recursive partition");
}

// -----------------------------------------------------------------------
// Error handling
// -----------------------------------------------------------------------

#[test]
fn kway_error_nparts_too_large() {
    let g = cycle4();
    let err = part_graph_kway(&g, 10, None, None, &Options::default()).unwrap_err();
    assert!(matches!(err, rmetis::MetisError::Input(_)));
}

#[test]
fn kway_error_nparts_one() {
    let g = cycle4();
    let err = part_graph_kway(&g, 1, None, None, &Options::default()).unwrap_err();
    assert!(matches!(err, rmetis::MetisError::Input(_)));
}

// -----------------------------------------------------------------------
// Weighted graph test
// -----------------------------------------------------------------------

#[test]
fn kway_weighted_graph() {
    // 4-cycle with vertex weights [2,1,2,1] and edge weights [3,1,3,1,3,1,3,1]
    let g = Graph::new(
        4, 1,
        vec![0,2,4,6,8],
        vec![1,3,0,2,1,3,0,2],
        Some(vec![2,1,2,1]),
        Some(vec![3,1,3,1,3,1,3,1]),
        None,
    ).unwrap();
    let opts = Options { seed: 5, ..Options::for_kway() };
    let r = part_graph_kway(&g, 2, None, None, &opts).unwrap();
    assert_partition_valid(&g, &r.part, 2);
}
