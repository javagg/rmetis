//! Graph data structure and operations (CSR format).

pub mod validate;

use crate::types::{Idx, MetisError, Real};
use validate::validate_csr;

/// Compressed Sparse Row graph representation.
///
/// Undirected graph: each edge (u, v) is stored twice — once in u's adjacency
/// list and once in v's adjacency list.
#[derive(Debug, Clone)]
pub struct Graph {
    /// Number of vertices
    pub nvtxs: usize,
    /// Number of undirected edges (adjncy.len() / 2)
    pub nedges: usize,
    /// CSR index array, length nvtxs + 1.
    /// Vertex v's neighbors are adjncy[xadj[v]..xadj[v+1]]
    pub xadj: Vec<Idx>,
    /// Adjacency array, length 2 * nedges
    pub adjncy: Vec<Idx>,
    /// Number of balancing constraints (≥ 1)
    pub ncon: usize,
    /// Vertex weights, length nvtxs * ncon (row-major: vwgt[v*ncon+c])
    /// None = uniform weight 1 per constraint
    pub vwgt: Option<Vec<Idx>>,
    /// Edge weights, length 2 * nedges (parallel to adjncy)
    /// None = uniform weight 1
    pub adjwgt: Option<Vec<Idx>>,
    /// Vertex communication sizes (for volume minimization)
    /// None = uniform size 1
    pub vsize: Option<Vec<Idx>>,
}

impl Graph {
    /// Construct and validate a graph from raw CSR arrays.
    pub fn new(
        nvtxs: usize,
        ncon: usize,
        xadj: Vec<Idx>,
        adjncy: Vec<Idx>,
        vwgt: Option<Vec<Idx>>,
        adjwgt: Option<Vec<Idx>>,
        vsize: Option<Vec<Idx>>,
    ) -> Result<Self, MetisError> {
        validate_csr(nvtxs, ncon, &xadj, &adjncy, vwgt.as_deref(), adjwgt.as_deref())?;
        let nedges = adjncy.len() / 2;
        Ok(Graph { nvtxs, nedges, ncon, xadj, adjncy, vwgt, adjwgt, vsize })
    }

    /// Create an unweighted graph (all weights = 1)
    pub fn new_unweighted(
        nvtxs: usize,
        xadj: Vec<Idx>,
        adjncy: Vec<Idx>,
    ) -> Result<Self, MetisError> {
        Self::new(nvtxs, 1, xadj, adjncy, None, None, None)
    }

    // -----------------------------------------------------------------------
    // Vertex accessors
    // -----------------------------------------------------------------------

    /// Neighbor vertex indices for vertex v
    #[inline]
    pub fn neighbors(&self, v: usize) -> &[Idx] {
        let s = self.xadj[v] as usize;
        let e = self.xadj[v + 1] as usize;
        &self.adjncy[s..e]
    }

    /// Edge weights for vertex v's adjacency list (parallel to neighbors)
    #[inline]
    pub fn edge_weights(&self, v: usize) -> Option<&[Idx]> {
        self.adjwgt.as_ref().map(|w| {
            let s = self.xadj[v] as usize;
            let e = self.xadj[v + 1] as usize;
            &w[s..e]
        })
    }

    /// Degree of vertex v
    #[inline]
    pub fn degree(&self, v: usize) -> usize {
        (self.xadj[v + 1] - self.xadj[v]) as usize
    }

    /// Weight of vertex v for constraint c (defaults to 1 if no vwgt)
    #[inline]
    pub fn vwgt_at(&self, v: usize, c: usize) -> Idx {
        match &self.vwgt {
            Some(w) => w[v * self.ncon + c],
            None => 1,
        }
    }

    /// Weight slice for vertex v (all constraints), defaults to &[1] if ncon=1
    pub fn vertex_weight_slice<'a>(&'a self, v: usize, scratch: &'a mut [Idx]) -> &'a [Idx] {
        match &self.vwgt {
            Some(w) => &w[v * self.ncon..(v + 1) * self.ncon],
            None => {
                scratch[..self.ncon].fill(1);
                &scratch[..self.ncon]
            }
        }
    }

    /// Edge weight between vertex v and its j-th neighbor (1-indexed into adjncy)
    #[inline]
    pub fn edge_weight_at(&self, pos: usize) -> Idx {
        match &self.adjwgt {
            Some(w) => w[pos],
            None => 1,
        }
    }

    // -----------------------------------------------------------------------
    // Partition-related queries
    // -----------------------------------------------------------------------

    /// Compute the edge cut of a given partition assignment.
    pub fn edge_cut(&self, part: &[Idx]) -> Idx {
        let mut cut = 0i32;
        for v in 0..self.nvtxs {
            let pv = part[v];
            let s = self.xadj[v] as usize;
            let e = self.xadj[v + 1] as usize;
            for j in s..e {
                let u = self.adjncy[j] as usize;
                if part[u] != pv {
                    cut += self.edge_weight_at(j);
                }
            }
        }
        cut / 2 // each edge counted twice
    }

    /// Compute per-partition weight vectors (shape: [nparts][ncon]).
    pub fn partition_weights(&self, part: &[Idx], nparts: usize) -> Vec<Vec<Idx>> {
        let mut pw = vec![vec![0i32; self.ncon]; nparts];
        for v in 0..self.nvtxs {
            let p = part[v] as usize;
            for c in 0..self.ncon {
                pw[p][c] += self.vwgt_at(v, c);
            }
        }
        pw
    }

    /// Total weight per constraint (length ncon).
    pub fn total_weight(&self) -> Vec<Idx> {
        let mut tw = vec![0i32; self.ncon];
        for v in 0..self.nvtxs {
            for c in 0..self.ncon {
                tw[c] += self.vwgt_at(v, c);
            }
        }
        tw
    }

    /// Collect boundary vertices: vertices with at least one neighbor in a
    /// different partition.
    pub fn boundary_vertices(&self, part: &[Idx]) -> Vec<usize> {
        (0..self.nvtxs)
            .filter(|&v| {
                let pv = part[v];
                self.neighbors(v).iter().any(|&u| part[u as usize] != pv)
            })
            .collect()
    }

    /// Check whether each partition is contiguous (connected subgraph).
    pub fn check_contiguous(&self, part: &[Idx], nparts: usize) -> bool {
        for p in 0..nparts as Idx {
            let verts: Vec<usize> = (0..self.nvtxs).filter(|&v| part[v] == p).collect();
            if verts.is_empty() {
                continue;
            }
            // BFS within partition p
            let mut visited = vec![false; self.nvtxs];
            let mut queue = std::collections::VecDeque::new();
            visited[verts[0]] = true;
            queue.push_back(verts[0]);
            let mut count = 1usize;
            while let Some(v) = queue.pop_front() {
                for &u in self.neighbors(v) {
                    let u = u as usize;
                    if !visited[u] && part[u] == p {
                        visited[u] = true;
                        count += 1;
                        queue.push_back(u);
                    }
                }
            }
            if count != verts.len() {
                return false;
            }
        }
        true
    }

    // -----------------------------------------------------------------------
    // Subgraph extraction (used in recursive bisection)
    // -----------------------------------------------------------------------

    /// Extract the subgraph induced by vertices with `part[v] == target_part`.
    /// Returns (subgraph, old_to_new vertex mapping, new_to_old vertex list).
    pub fn extract_subgraph(
        &self,
        part: &[Idx],
        target_part: Idx,
    ) -> (Graph, Vec<usize>, Vec<usize>) {
        let new_to_old: Vec<usize> = (0..self.nvtxs)
            .filter(|&v| part[v] == target_part)
            .collect();
        let snvtxs = new_to_old.len();

        let mut old_to_new = vec![usize::MAX; self.nvtxs];
        for (new, &old) in new_to_old.iter().enumerate() {
            old_to_new[old] = new;
        }

        let mut xadj = vec![0i32; snvtxs + 1];
        for (new, &old) in new_to_old.iter().enumerate() {
            let deg = self.neighbors(old)
                .iter()
                .filter(|&&u| part[u as usize] == target_part)
                .count();
            xadj[new + 1] = xadj[new] + deg as Idx;
        }

        let nedges_stored = xadj[snvtxs] as usize;
        let mut adjncy = Vec::with_capacity(nedges_stored);
        let mut adjwgt_vec: Option<Vec<Idx>> = self.adjwgt.as_ref().map(|_| Vec::with_capacity(nedges_stored));

        for (new, &old) in new_to_old.iter().enumerate() {
            let s = self.xadj[old] as usize;
            let e = self.xadj[old + 1] as usize;
            for j in s..e {
                let u = self.adjncy[j] as usize;
                if part[u] == target_part {
                    adjncy.push(old_to_new[u] as Idx);
                    if let Some(ref mut aw) = adjwgt_vec {
                        aw.push(self.edge_weight_at(j));
                    }
                }
            }
            let _ = new; // suppress warning
        }

        // Extract vertex weights
        let vwgt_sub = self.vwgt.as_ref().map(|vw| {
            new_to_old.iter().flat_map(|&old| {
                vw[old * self.ncon..(old + 1) * self.ncon].iter().copied()
            }).collect::<Vec<_>>()
        });

        let vsize_sub = self.vsize.as_ref().map(|vs| {
            new_to_old.iter().map(|&old| vs[old]).collect::<Vec<_>>()
        });

        let g = Graph {
            nvtxs: snvtxs,
            nedges: adjncy.len() / 2,
            ncon: self.ncon,
            xadj,
            adjncy,
            vwgt: vwgt_sub,
            adjwgt: adjwgt_vec,
            vsize: vsize_sub,
        };
        (g, old_to_new, new_to_old)
    }

    // -----------------------------------------------------------------------
    // Coarse-graph construction (used during multilevel coarsening)
    // -----------------------------------------------------------------------

    /// Build a coarse graph from a matching.
    ///
    /// `matching[v]` = vertex matched with v (or `v` itself if unmatched).
    /// Returns (coarse_graph, cmap) where cmap[v] = super-vertex id.
    pub fn build_coarse_graph(&self, matching: &[usize]) -> (Graph, Vec<usize>) {
        // Assign super-vertex IDs
        let mut cmap = vec![usize::MAX; self.nvtxs];
        let mut cnvtxs = 0usize;
        for v in 0..self.nvtxs {
            if cmap[v] == usize::MAX {
                cmap[v] = cnvtxs;
                let u = matching[v];
                if u != v && cmap[u] == usize::MAX {
                    cmap[u] = cnvtxs;
                }
                cnvtxs += 1;
            }
        }

        // Aggregate vertex weights
        let mut cvwgt: Option<Vec<Idx>> = self.vwgt.as_ref().map(|_| vec![0i32; cnvtxs * self.ncon]);
        let mut cvsize: Option<Vec<Idx>> = self.vsize.as_ref().map(|_| vec![0i32; cnvtxs]);

        for v in 0..self.nvtxs {
            let cv = cmap[v];
            for c in 0..self.ncon {
                if let Some(ref mut cw) = cvwgt {
                    cw[cv * self.ncon + c] += self.vwgt_at(v, c);
                }
            }
            if let Some(ref mut cs) = cvsize {
                let vs = self.vsize.as_ref().map(|s| s[v]).unwrap_or(1);
                cs[cv] += vs;
            }
        }

        // Build coarse adjacency (merge multi-edges by summing weights)
        // Use a temporary dense row buffer per coarse vertex
        let mut cxadj = vec![0i32; cnvtxs + 1];
        // First pass: collect edges into per-super-vertex edge lists
        // We use a Vec<Vec<(usize, Idx)>> for simplicity at this stage
        let mut cadj_lists: Vec<Vec<(usize, Idx)>> = vec![Vec::new(); cnvtxs];

        for v in 0..self.nvtxs {
            let cv = cmap[v];
            let s = self.xadj[v] as usize;
            let e = self.xadj[v + 1] as usize;
            for j in s..e {
                let u = self.adjncy[j] as usize;
                let cu = cmap[u];
                if cu == cv {
                    continue; // internal edge, skip
                }
                let ew = self.edge_weight_at(j);
                // Merge: find existing edge to cu
                if let Some(existing) = cadj_lists[cv].iter_mut().find(|(n, _)| *n == cu) {
                    existing.1 += ew;
                } else {
                    cadj_lists[cv].push((cu, ew));
                }
            }
        }

        // Build CSR from edge lists
        for cv in 0..cnvtxs {
            cxadj[cv + 1] = cxadj[cv] + cadj_lists[cv].len() as Idx;
        }
        let cnedges_stored = cxadj[cnvtxs] as usize;
        let mut cadjncy = Vec::with_capacity(cnedges_stored);
        let mut cadjwgt: Option<Vec<Idx>> = if self.adjwgt.is_some() {
            Some(Vec::with_capacity(cnedges_stored))
        } else {
            // still track weights if any edge had weight
            None
        };
        // Always track edge weights in coarse graph for refinement quality
        let track_weights = true;
        let mut cadjwgt_vec = Vec::with_capacity(cnedges_stored);

        for cv in 0..cnvtxs {
            for (cu, ew) in &cadj_lists[cv] {
                cadjncy.push(*cu as Idx);
                cadjwgt_vec.push(*ew);
            }
        }

        let has_weights = self.adjwgt.is_some()
            || cadjwgt_vec.iter().any(|&w| w != 1);
        let _ = (cadjwgt, track_weights);
        let cadjwgt_out = if has_weights { Some(cadjwgt_vec) } else { None };

        let cnedges = cadjncy.len() / 2;

        let cg = Graph {
            nvtxs: cnvtxs,
            nedges: cnedges,
            ncon: self.ncon,
            xadj: cxadj,
            adjncy: cadjncy,
            vwgt: cvwgt,
            adjwgt: cadjwgt_out,
            vsize: cvsize,
        };
        (cg, cmap)
    }

    // -----------------------------------------------------------------------
    // Target weight utilities
    // -----------------------------------------------------------------------

    /// Compute target weights for each partition given `tpwgts`.
    /// Returns a Vec of shape [nparts][ncon] scaled to integer vertex weights.
    pub fn compute_target_weights(
        &self,
        nparts: usize,
        tpwgts: Option<&[Real]>,
    ) -> Vec<Vec<f64>> {
        let total = self.total_weight();
        let mut targets = vec![vec![0.0f64; self.ncon]; nparts];
        for p in 0..nparts {
            for c in 0..self.ncon {
                let t = if let Some(tp) = tpwgts {
                    tp[p * self.ncon + c] as f64
                } else {
                    1.0 / nparts as f64
                };
                targets[p][c] = t * total[c] as f64;
            }
        }
        targets
    }
}

impl core::fmt::Display for Graph {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "Graph(nvtxs={}, nedges={}, ncon={})",
            self.nvtxs, self.nedges, self.ncon
        )
    }
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// 4-vertex cycle: 0-1-2-3-0
    fn cycle4() -> Graph {
        Graph::new_unweighted(
            4,
            vec![0, 2, 4, 6, 8],
            vec![1, 3, 0, 2, 1, 3, 0, 2],
        ).unwrap()
    }

    #[test]
    fn basic_construction() {
        let g = cycle4();
        assert_eq!(g.nvtxs, 4);
        assert_eq!(g.nedges, 4);
    }

    #[test]
    fn neighbors() {
        let g = cycle4();
        assert_eq!(g.neighbors(0), &[1, 3]);
        assert_eq!(g.neighbors(1), &[0, 2]);
    }

    #[test]
    fn edge_cut_partition() {
        let g = cycle4();
        // partition: {0,1} | {2,3}
        let part = vec![0, 0, 1, 1];
        // Edges crossing: 1-2, 3-0 → 2 edges
        assert_eq!(g.edge_cut(&part), 2);
    }

    #[test]
    fn partition_weights_uniform() {
        let g = cycle4();
        let part = vec![0, 0, 1, 1];
        let pw = g.partition_weights(&part, 2);
        assert_eq!(pw[0], vec![2]);
        assert_eq!(pw[1], vec![2]);
    }

    #[test]
    fn total_weight_uniform() {
        let g = cycle4();
        assert_eq!(g.total_weight(), vec![4]);
    }

    #[test]
    fn boundary_vertices() {
        let g = cycle4();
        let part = vec![0, 0, 1, 1];
        let mut bv = g.boundary_vertices(&part);
        bv.sort();
        assert_eq!(bv, vec![0, 1, 2, 3]); // all 4 are on boundary
    }

    #[test]
    fn extract_subgraph() {
        let g = cycle4();
        let part = vec![0, 0, 1, 1];
        let (sg, _o2n, n2o) = g.extract_subgraph(&part, 0);
        assert_eq!(sg.nvtxs, 2);
        assert_eq!(n2o, vec![0, 1]);
        // Only internal edge 0-1 remains
        assert_eq!(sg.nedges, 1);
    }

    #[test]
    fn build_coarse_graph_path() {
        // Path: 0-1-2-3
        let g = Graph::new_unweighted(
            4,
            vec![0, 1, 3, 5, 6],
            vec![1, 0, 2, 1, 3, 2],
        ).unwrap();
        // Match 0-1, 2-3
        let matching = vec![1, 0, 3, 2];
        let (cg, cmap) = g.build_coarse_graph(&matching);
        assert_eq!(cg.nvtxs, 2);
        assert_eq!(cg.nedges, 1);
        assert_eq!(cmap[0], cmap[1]);
        assert_eq!(cmap[2], cmap[3]);
        assert_ne!(cmap[0], cmap[2]);
    }
}
