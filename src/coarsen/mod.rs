//! Multilevel coarsening framework.

pub mod rm;
pub mod shem;

use crate::comm::{Comm, SingleComm};
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
    let single = SingleComm;
    coarsen_with_comm(graph, nparts, options, rng, &single)
}

/// Coarsen with optional MPI parallelism.
///
/// In multi-rank mode, each rank computes matching proposals only for its
/// owned vertices (`v % nprocs == rank`). Proposals are gathered via
/// `all_gather` and conflicts are resolved deterministically (the lowest
/// rank that proposed a match for a vertex wins). The resulting matching
/// is identical on all ranks, so `build_coarse_graph` is run locally.
pub fn coarsen_with_comm(
    graph: &Graph,
    nparts: usize,
    options: &Options,
    rng: &mut impl Rng,
    comm: &dyn Comm,
) -> Vec<CoarseLevel> {
    let coarsen_limit = (20 * nparts).max(20);
    let nprocs = comm.size().max(1) as usize;
    let rank = comm.rank() as usize;
    let mut levels: Vec<CoarseLevel> = Vec::new();
    let mut current = graph.clone();

    loop {
        let matching = if nprocs <= 1 {
            // Serial path — same as before
            match options.ctype {
                CType::Rm   => rm::random_matching(&current, rng),
                CType::Shem => shem::shem_matching(&current, rng),
            }
        } else {
            parallel_matching(&current, options, rng, rank, nprocs, comm)
        };

        let (cg, cmap) = current.build_coarse_graph(&matching);
        let ratio = cg.nvtxs as f32 / current.nvtxs as f32;

        levels.push(CoarseLevel { graph: cg.clone(), cmap });

        if cg.nvtxs <= coarsen_limit || ratio > 0.95 {
            break;
        }
        current = cg;
    }
    levels
}

/// Distributed matching: each rank proposes matches for its owned vertices.
/// Proposals are exchanged via `all_gather`; conflicts resolved by lowest rank.
fn parallel_matching(
    graph: &Graph,
    options: &Options,
    rng: &mut impl Rng,
    rank: usize,
    nprocs: usize,
    comm: &dyn Comm,
) -> Vec<usize> {
    let n = graph.nvtxs;

    // Start with full local matching for the entire graph (each rank does this
    // independently with a rank-offset seed to get diverse proposals).
    let full_matching = match options.ctype {
        CType::Rm   => rm::random_matching(graph, rng),
        CType::Shem => shem::shem_matching(graph, rng),
    };

    // Each rank keeps only proposals for vertices it owns
    // Encode as pairs [v, partner, v, partner, ...]
    let local_proposals: Vec<i32> = (0..n)
        .filter(|&v| v % nprocs == rank)
        .flat_map(|v| [v as i32, full_matching[v] as i32])
        .collect();

    // Exchange proposals
    let all_proposals = comm.all_gather_i32_vec(&local_proposals);

    // Build a "claimed_by" array: for each vertex, which rank first claimed it
    // (to resolve conflicts). Final matching is self-match by default.
    let mut matching: Vec<usize> = (0..n).collect();
    // For each vertex, track which rank's proposal we're accepting
    let mut claimed_by: Vec<Option<usize>> = vec![None; n];
    let mut claimed_partner: Vec<usize> = (0..n).collect(); // the proposed partner

    for (proposing_rank, proposals) in all_proposals.iter().enumerate() {
        let mut i = 0;
        while i + 1 < proposals.len() {
            let v = proposals[i] as usize;
            let u = proposals[i + 1] as usize;
            i += 2;

            if v == u {
                // self-match, skip
                continue;
            }

            // Accept this proposal only if neither v nor u has been claimed
            let v_free = claimed_by[v].is_none();
            let u_free = claimed_by[u].is_none();

            if v_free && u_free {
                claimed_by[v] = Some(proposing_rank);
                claimed_by[u] = Some(proposing_rank);
                claimed_partner[v] = u;
                claimed_partner[u] = v;
            }
        }
    }

    // Finalise matching: use claimed_partner where a match was accepted
    for v in 0..n {
        if claimed_by[v].is_some() {
            matching[v] = claimed_partner[v];
        }
        // else: matching[v] == v (self-match, already set)
    }

    matching
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::comm::SingleComm;
    use crate::graph::Graph;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn path6() -> Graph {
        // 0-1-2-3-4-5 undirected path
        Graph::new_unweighted(
            6,
            vec![0, 1, 3, 5, 7, 9, 10],
            vec![1, 0, 2, 1, 3, 2, 4, 3, 5, 4],
        ).unwrap()
    }

    #[test]
    fn coarsen_with_comm_single_rank_matches_serial() {
        let g = path6();
        let opts = Options::default();
        let mut rng_a = SmallRng::seed_from_u64(99);
        let mut rng_b = SmallRng::seed_from_u64(99);
        let levels_a = coarsen(&g, 2, &opts, &mut rng_a);
        let levels_b = coarsen_with_comm(&g, 2, &opts, &mut rng_b, &SingleComm);
        assert_eq!(levels_a.len(), levels_b.len(), "coarsen levels differ");
        for (la, lb) in levels_a.iter().zip(levels_b.iter()) {
            assert_eq!(la.graph.nvtxs, lb.graph.nvtxs);
        }
    }
}

