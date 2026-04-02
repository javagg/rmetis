//! Initial partitioning algorithms.

pub mod grow;
pub mod random;

use crate::graph::Graph;
use crate::types::{Idx, IpType, Options, Real};
use rand::Rng;

/// Produce an initial partition of `graph` into `nparts` parts.
pub fn initial_partition(
    graph: &Graph,
    nparts: usize,
    tpwgts: &[f64],
    options: &Options,
    rng: &mut impl Rng,
) -> Vec<Idx> {
    match options.iptype {
        IpType::Grow | IpType::MetisRb => grow::grow_partition(graph, nparts, tpwgts, rng),
        IpType::Random => random::random_partition(graph, nparts, rng),
    }
}
