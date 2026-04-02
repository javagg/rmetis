//! WebAssembly bindings via wasm-bindgen.

use wasm_bindgen::prelude::*;
use crate::graph::Graph;
use crate::partition::{kway, recursive, nd};
use crate::types::{Idx, Options, Real};

/// A graph object exposed to JavaScript.
#[wasm_bindgen]
pub struct WasmGraph {
    inner: Graph,
}

#[wasm_bindgen]
impl WasmGraph {
    /// Create a new unweighted graph.
    ///
    /// `xadj` has length `nvtxs + 1`, `adjncy` has length `xadj[nvtxs]`.
    #[wasm_bindgen(constructor)]
    pub fn new(
        nvtxs: usize,
        xadj: Vec<i32>,
        adjncy: Vec<i32>,
    ) -> Result<WasmGraph, JsValue> {
        Graph::new_unweighted(nvtxs, xadj, adjncy)
            .map(|g| WasmGraph { inner: g })
            .map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Set vertex weights (length = nvtxs * ncon).
    pub fn set_vertex_weights(&mut self, vwgt: Vec<i32>, ncon: usize) -> Result<(), JsValue> {
        if vwgt.len() != self.inner.nvtxs * ncon {
            return Err(JsValue::from_str("vwgt length mismatch"));
        }
        self.inner.vwgt = Some(vwgt);
        self.inner.ncon = ncon;
        Ok(())
    }

    /// Set edge weights (length = adjncy.len()).
    pub fn set_edge_weights(&mut self, adjwgt: Vec<i32>) -> Result<(), JsValue> {
        if adjwgt.len() != self.inner.adjncy.len() {
            return Err(JsValue::from_str("adjwgt length mismatch"));
        }
        self.inner.adjwgt = Some(adjwgt);
        Ok(())
    }

    /// Number of vertices.
    pub fn nvtxs(&self) -> usize { self.inner.nvtxs }

    /// Number of undirected edges.
    pub fn nedges(&self) -> usize { self.inner.nedges }
}

/// Options object for partitioning.
#[wasm_bindgen]
pub struct WasmOptions {
    inner: Options,
}

#[wasm_bindgen]
impl WasmOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> WasmOptions {
        WasmOptions { inner: Options::default() }
    }

    pub fn set_seed(&mut self, seed: u32) { self.inner.seed = seed as u64; }
    pub fn set_ncuts(&mut self, n: usize) { self.inner.ncuts = n; }
    pub fn set_niter(&mut self, n: usize) { self.inner.niter = n; }
    pub fn set_ufactor(&mut self, u: u32) { self.inner.ufactor = u; }
    pub fn set_contig(&mut self, v: bool) { self.inner.contig = v; }
    pub fn set_minconn(&mut self, v: bool) { self.inner.minconn = v; }
}

/// Result of a graph partitioning.
#[wasm_bindgen]
pub struct WasmPartitionResult {
    part: Vec<i32>,
    objval: i32,
}

#[wasm_bindgen]
impl WasmPartitionResult {
    pub fn part(&self) -> Vec<i32> { self.part.clone() }
    pub fn objval(&self) -> i32 { self.objval }
}

/// Result of nested dissection.
#[wasm_bindgen]
pub struct WasmNDResult {
    perm: Vec<i32>,
    iperm: Vec<i32>,
}

#[wasm_bindgen]
impl WasmNDResult {
    pub fn perm(&self) -> Vec<i32> { self.perm.clone() }
    pub fn iperm(&self) -> Vec<i32> { self.iperm.clone() }
}

/// K-way graph partitioning.
#[wasm_bindgen]
pub fn part_graph_kway(
    graph: &WasmGraph,
    nparts: usize,
    options: Option<WasmOptions>,
) -> Result<WasmPartitionResult, JsValue> {
    let opts = options.map(|o| o.inner).unwrap_or_default();
    kway::partition_kway(&graph.inner, nparts, None, None, &opts)
        .map(|r| WasmPartitionResult { part: r.part, objval: r.objval })
        .map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Recursive bisection graph partitioning.
#[wasm_bindgen]
pub fn part_graph_recursive(
    graph: &WasmGraph,
    nparts: usize,
    options: Option<WasmOptions>,
) -> Result<WasmPartitionResult, JsValue> {
    let opts = options.map(|o| o.inner).unwrap_or_default();
    recursive::partition_recursive(&graph.inner, nparts, None, None, &opts)
        .map(|r| WasmPartitionResult { part: r.part, objval: r.objval })
        .map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Nested dissection ordering.
#[wasm_bindgen]
pub fn node_nd(
    graph: &WasmGraph,
    options: Option<WasmOptions>,
) -> Result<WasmNDResult, JsValue> {
    let opts = options.map(|o| o.inner).unwrap_or_default();
    nd::node_nd(&graph.inner, &opts)
        .map(|r| WasmNDResult { perm: r.perm, iperm: r.iperm })
        .map_err(|e| JsValue::from_str(&e.to_string()))
}
