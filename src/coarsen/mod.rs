//! Multilevel coarsening framework.

pub mod rm;
pub mod shem;

use crate::graph::Graph;
use crate::types::{CType, Options};
use rand::Rng;

/// A single level in the coarsening hierarchy.
pub struct CoarseLevel {
    /// The coarsened graph at this level
    pub graph: Graph,
    /// Maps fine-level vertex → super-vertex in this coarse graph
    pub cmap: Vec<usize>,
}

/// Coarsen `graph` until it is small enough for direct partitioning.
///
/// Returns levels from finest-1 to coarsest (level 0 is the coarsened
/// version of the original; the last entry is the coarsest graph).
/// The original graph is NOT included.
pub fn coarsen(graph: &Graph, nparts: usize, options: &Options, rng: &mut impl Rng) -> Vec<CoarseLevel> {
    let coarsen_limit = (20 * nparts).max(20);
    let mut levels: Vec<CoarseLevel> = Vec::new();
    let mut current = graph.clone();

    loop {
        let matching = match options.ctype {
            CType::Rm   => rm::random_matching(&current, rng),
            CType::Shem => shem::shem_matching(&current, rng),
        };

        let (cg, cmap) = current.build_coarse_graph(&matching);
        let ratio = cg.nvtxs as f32 / current.nvtxs as f32;

        levels.push(CoarseLevel { graph: cg.clone(), cmap });

        // Stop conditions
        if cg.nvtxs <= coarsen_limit || ratio > 0.95 {
            break;
        }
        current = cg;
    }
    levels
}
