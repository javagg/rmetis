//! Refinement algorithms.

pub mod greedy;
pub mod fm;
pub mod fm_kway;

use crate::graph::Graph;
use crate::types::{Idx, Options, Real, RType};

/// Refine a 2-way partition in-place; returns new edge cut.
pub fn refine_2way(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],
    ubvec: &[f64],
    options: &Options,
) -> Idx {
    match options.rtype {
        RType::Greedy => greedy::greedy_refine(graph, part, nparts, tpwgts, ubvec, options.niter),
        RType::Fm     => fm::fm_refine(graph, part, tpwgts, ubvec, options.niter),
    }
}

/// Refine a k-way partition in-place; returns new edge cut.
pub fn refine_kway(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],
    ubvec: &[f64],
    options: &Options,
) -> Idx {
    match options.rtype {
        RType::Greedy => greedy::greedy_refine(graph, part, nparts, tpwgts, ubvec, options.niter),
        RType::Fm     => fm_kway::fm_kway_refine(graph, part, nparts, tpwgts, ubvec, options.niter),
    }
}
