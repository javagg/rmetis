use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rmetis::{Graph, Options, part_graph_kway, part_graph_recursive};

fn build_grid(rows: usize, cols: usize) -> Graph {
    let n = rows * cols;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if c + 1 < cols { adj[v].push(v + 1); adj[v + 1].push(v); }
            if r + 1 < rows { adj[v].push(v + cols); adj[v + cols].push(v); }
        }
    }
    let mut xadj = vec![0i32; n + 1];
    let mut adjncy = Vec::new();
    for v in 0..n {
        adj[v].sort();
        adj[v].dedup();
        xadj[v + 1] = xadj[v] + adj[v].len() as i32;
        for &u in &adj[v] { adjncy.push(u as i32); }
    }
    Graph::new_unweighted(n, xadj, adjncy).unwrap()
}

fn bench_kway(c: &mut Criterion) {
    let graphs = [
        ("grid10x10", build_grid(10, 10)),
        ("grid20x20", build_grid(20, 20)),
        ("grid50x50", build_grid(50, 50)),
    ];

    let mut group = c.benchmark_group("PartGraphKway");
    for (name, graph) in &graphs {
        for nparts in [2, 4, 8] {
            let opts = Options { seed: 42, ..Options::for_kway() };
            group.bench_with_input(
                BenchmarkId::new(*name, nparts),
                &nparts,
                |b, &k| b.iter(|| part_graph_kway(graph, k, None, None, &opts).unwrap()),
            );
        }
    }
    group.finish();
}

fn bench_recursive(c: &mut Criterion) {
    let graph = build_grid(20, 20);
    let mut group = c.benchmark_group("PartGraphRecursive");
    for nparts in [2, 4, 8] {
        let opts = Options { seed: 42, ..Options::for_recursive() };
        group.bench_with_input(
            BenchmarkId::from_parameter(nparts),
            &nparts,
            |b, &k| b.iter(|| part_graph_recursive(&graph, k, None, None, &opts).unwrap()),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_kway, bench_recursive);
criterion_main!(benches);
