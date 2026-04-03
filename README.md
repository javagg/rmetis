# rmetis

Pure Rust graph partitioning library, targeting [METIS 5.1.x](https://karypis.github.io/glaros/software/metis/overview.html) API compatibility. Compiles to WebAssembly with no native dependencies. Supports MPI-level distributed parallelism via [jsmpi](https://github.com/jsmpi/jsmpi) (WASM) and the [mpi](https://crates.io/crates/mpi) crate (native).

## Features

- **Pure Rust** — no C/C++ dependencies, no unsafe code outside the optional C ABI layer
- **WASM-ready** — targets `wasm32-unknown-unknown` out of the box
- **METIS-compatible** — same three core functions and option semantics as METIS 5.1.x
- **Deterministic** — fixed seed produces identical results across platforms
- **Multi-constraint** — supports balancing multiple weight dimensions simultaneously
- **MPI parallel** — ParMETIS-style distributed parallelism: parallel ncuts, distributed coarsening, distributed greedy refinement

## Algorithms

| Phase | Options |
|-------|---------|
| Coarsening | Random Matching (RM), Sorted Heavy-Edge Matching (SHEM) |
| Initial partitioning | Graph-growing (BFS), Random |
| Refinement | Greedy boundary (with balance-repair), Fiduccia-Mattheyses (FM), K-way FM |

## Quick Start

```toml
[dependencies]
rmetis = "0.1"
```

```rust
use rmetis::{Graph, Options, part_graph_kway};

// Build a 4-vertex cycle: 0-1-2-3-0
let graph = Graph::new_unweighted(
    4,
    vec![0, 2, 4, 6, 8],           // xadj
    vec![1, 3, 0, 2, 1, 3, 0, 2],  // adjncy (each undirected edge stored twice)
)?;

let result = part_graph_kway(&graph, 2, None, None, &Options::default())?;
println!("partition: {:?}", result.part);   // e.g. [0, 0, 1, 1]
println!("edge cut:  {}", result.objval);   // 2
```

## API

### Graph construction

```rust
// Unweighted graph
let g = Graph::new_unweighted(nvtxs, xadj, adjncy)?;

// Weighted graph (vertex weights + edge weights)
let g = Graph::new(nvtxs, ncon, xadj, adjncy, Some(vwgt), Some(adjwgt), None)?;
```

Graphs use **CSR (Compressed Sparse Row)** format — the same layout as METIS:

```
xadj[v]..xadj[v+1]  →  range in adjncy for vertex v's neighbors
```

Undirected graphs store each edge twice (once per direction).

### Partitioning

```rust
// K-way partitioning (faster, good for large k)
let result = part_graph_kway(&graph, nparts, tpwgts, ubvec, &options)?;

// Recursive bisection (higher quality, better for small k)
let result = part_graph_recursive(&graph, nparts, tpwgts, ubvec, &options)?;

// result.part[v]   — partition ID for vertex v (0..nparts)
// result.objval    — edge cut
```

### Nested dissection ordering

```rust
// Fill-reducing vertex ordering for sparse matrix factorization
let result = node_nd(&graph, &options)?;

// result.perm[v]    — new position of original vertex v
// result.iperm[i]   — original vertex at new position i
```

### Options

```rust
let options = Options {
    seed:     42,           // 0 = entropy-seeded
    ncuts:    3,            // try N times, keep best result
    niter:    10,           // refinement iterations per level
    ufactor:  30,           // imbalance tolerance × 1000 (30 = 3%)
    ctype:    CType::Shem,  // coarsening: Shem | Rm
    rtype:    RType::Fm,    // refinement: Fm | Greedy
    ..Default::default()
};
```

### MPI parallel partitioning

```rust
use rmetis::comm;
use rmetis::partition::kway::partition_kway;

// Initialise the MPI world communicator (auto-selects backend)
let world = comm::world();  // Box<dyn Comm>

// With ncuts=8 and 4 MPI ranks, each rank handles 2 cut attempts.
// The globally best partition is broadcast back to every rank.
let options = Options { seed: 1, ncuts: 8, ..Options::for_kway() };
let result = partition_kway(&graph, 4, None, None, &options, Some(world.as_ref()))?;
```

`comm::world()` returns a `SingleComm` (no-op) when only one process is present, so the same code works both serially and in parallel.

## Cargo Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | ✓ | Standard library support |
| `c-api` | | C ABI functions (`RMETIS_PartGraphKway`, etc.) |
| `wasm` | | WebAssembly bindings via `wasm-bindgen` + jsmpi MPI backend |
| `parallel` | | Shared-memory parallelism via `rayon` |

MPI support is automatic via target-conditional dependencies:

| Target | MPI backend |
|--------|-------------|
| `wasm32-*` | [jsmpi](https://github.com/jsmpi/jsmpi) (git submodule at `vendor/jsmpi`) |
| all others | [mpi](https://crates.io/crates/mpi) 0.8 (wraps OpenMPI / MPICH) |

## Parallel Algorithms

### Parallel ncuts

When `comm.size() > 1`, the `ncuts` independent partition attempts are distributed across MPI ranks. Rank `r` handles attempts `r, r+P, r+2P, …` with different seeds. All ranks `all_reduce_min` the best edge cut and `broadcast` the winning partition.

### Distributed coarsening

Each rank computes a matching for its owned vertices (`v % P == rank`) independently. Proposals are exchanged via `all_gather`; conflicts (two ranks claiming the same vertex) are resolved by lowest-rank-wins. All ranks end up with the identical matching and build the coarse graph locally.

### Distributed greedy refinement

Boundary vertices are split across ranks (`v % P == rank`). Each rank runs a **balance-repair pass** first (moves vertices out of overweight partitions even at neutral gain), then a **cut-optimisation pass** (positive-gain moves only). After each pass all ranks exchange proposed moves via `all_gather` and apply them atomically.

## WebAssembly

### Build

```bash
# Install wasm-pack
cargo install wasm-pack

# Build for browser
wasm-pack build --target web --release -- --features wasm

# Build for Node.js
wasm-pack build --target nodejs --release -- --features wasm
```

### JavaScript usage

```javascript
import init, { WasmGraph, part_graph_kway } from './pkg/rmetis.js';

await init();

const graph = new WasmGraph(
    4,
    new Int32Array([0, 2, 4, 6, 8]),
    new Int32Array([1, 3, 0, 2, 1, 3, 0, 2])
);

const result = part_graph_kway(graph, 2);
console.log(result.part());    // Int32Array
console.log(result.objval());  // number
```

## C ABI (METIS Drop-in)

Build as a shared library with `--features c-api`:

```bash
cargo build --release --features c-api
# → target/release/librmetis.so (Linux)
# → target/release/rmetis.dll   (Windows)
```

The exported symbols match METIS 5.1.x signatures with an `RMETIS_` prefix:

```c
#include "include/rmetis.h"

idx_t options[RMETIS_NOPTIONS];
RMETIS_SetDefaultOptions(options);

idx_t objval;
RMETIS_PartGraphKway(
    &nvtxs, &ncon, xadj, adjncy,
    NULL, NULL, NULL,           // vwgt, vsize, adjwgt (NULL = uniform)
    &nparts, NULL, NULL,        // tpwgts, ubvec (NULL = defaults)
    options, &objval, part
);
```

## Graph Format

Input graphs use 0-based CSR. To load METIS `.graph` files (1-based), subtract 1 from all indices.

```
# METIS .graph format
nvtxs nedges [fmt]
<neighbors of vertex 1, space-separated>
<neighbors of vertex 2>
...
```

## Performance

Indicative timings on a single core (x86_64, release build):

| Graph | Vertices | k | Algorithm | Time |
|-------|----------|---|-----------|------|
| 10×10 grid | 100 | 4 | K-way | ~1 ms |
| 20×20 grid | 400 | 8 | K-way | ~5 ms |
| 50×50 grid | 2500 | 16 | K-way | ~50 ms |

## Project Structure

```
src/
├── lib.rs              # Public API
├── types.rs            # Types: Idx, Real, Options, MetisError
├── comm.rs             # MPI comm abstraction (SingleComm / jsmpi / mpi)
├── graph/              # CSR graph structure and operations
├── coarsen/            # RM and SHEM coarsening (serial + distributed)
├── initial/            # BFS-grow and random initial partitioning
├── refine/             # Greedy (with balance-repair), FM, and K-way FM
├── partition/          # kway (parallel ncuts), recursive bisection, NodeND
├── separator/          # Vertex separator (used by NodeND)
└── ffi/                # C ABI and WASM bindings
vendor/
└── jsmpi/              # Git submodule — MPI compatibility layer for WASM
tests/
├── correctness.rs      # Algorithm correctness and balance tests
└── parallel.rs         # MPI abstraction and parallel algorithm tests
```

## References

1. Karypis & Kumar (1998). *A fast and high quality multilevel scheme for partitioning irregular graphs.* SIAM J. Sci. Comput., 20(1).
2. Karypis & Kumar (1998). *Multilevel k-way partitioning scheme for irregular graphs.* J. Parallel Distrib. Comput., 48(1).
3. Karypis, Schloegel & Kumar (2003). *ParMETIS: Parallel Graph Partitioning and Sparse Matrix Ordering Library.* Technical Report.
4. [METIS 5.1.0 Manual](https://karypis.github.io/glaros/files/sw/metis/manual.pdf)

## License

MIT
