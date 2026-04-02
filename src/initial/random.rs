//! Random initial partitioning.

use crate::graph::Graph;
use crate::types::Idx;
use rand::Rng;

/// Randomly assign vertices to partitions (uniform distribution).
pub fn random_partition(graph: &Graph, nparts: usize, rng: &mut impl Rng) -> Vec<Idx> {
    (0..graph.nvtxs).map(|_| rng.gen_range(0..nparts) as Idx).collect()
}
