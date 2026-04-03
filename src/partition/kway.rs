//! Multilevel k-way partitioning.

use crate::coarsen;
use crate::comm::{Comm, SingleComm};
use crate::graph::Graph;
use crate::initial;
use crate::refine;
use crate::refine::greedy::greedy_refine_with_comm;
use crate::types::{Idx, MetisError, Options, PartitionResult, Real, RType};
use rand::rngs::SmallRng;
use rand::SeedableRng;

/// Partition `graph` into `nparts` partitions using the multilevel k-way scheme.
///
/// When `comm` has `size() > 1`, the `ncuts` attempts are distributed across
/// MPI ranks (each rank handles roughly `ncuts / size` attempts with rank-offset
/// seeds), and the globally best partition is broadcast back to every rank.
/// Within each cut attempt the distributed coarsening and greedy refinement
/// kernels are also used (for FM refinement, falls back to serial).
pub fn partition_kway(
    graph: &Graph,
    nparts: usize,
    tpwgts: Option<&[Real]>,
    ubvec: Option<&[Real]>,
    options: &Options,
    comm: Option<&dyn Comm>,
) -> Result<PartitionResult, MetisError> {
    let single = SingleComm;
    let comm: &dyn Comm = comm.unwrap_or(&single);

    if nparts < 2 {
        return Err(MetisError::Input(format!("nparts must be >= 2, got {}", nparts)));
    }
    if nparts > graph.nvtxs {
        return Err(MetisError::Input(format!(
            "nparts={} > nvtxs={}", nparts, graph.nvtxs
        )));
    }

    let ncon = graph.ncon;
    let rank = comm.rank() as usize;
    let nprocs = comm.size().max(1) as usize;

    let ub: Vec<f64> = if let Some(ub) = ubvec {
        ub.iter().map(|&x| x as f64).collect()
    } else {
        let default_ub = 1.0 + options.ufactor as f64 / 1000.0;
        vec![default_ub; ncon]
    };

    let mut best_result: Option<PartitionResult> = None;

    // Each rank handles cut attempts with indices ≡ rank (mod nprocs).
    for cut_attempt in (rank..options.ncuts).step_by(nprocs) {
        let seed = if options.seed == 0 {
            let mut entropy_rng = SmallRng::from_entropy();
            use rand::Rng;
            entropy_rng.gen::<u64>()
        } else {
            options.seed.wrapping_add(cut_attempt as u64)
        };
        let mut rng = SmallRng::seed_from_u64(seed);

        let levels = coarsen::coarsen_with_comm(graph, nparts, options, &mut rng, &SingleComm);

        let coarsest = levels.last().map(|l| &l.graph).unwrap_or(graph);
        let coarse_targets = coarsest.compute_target_weights(nparts, tpwgts);
        let coarse_flat: Vec<f64> = coarse_targets.iter()
            .flat_map(|r| r.iter().copied()).collect();
        let mut coarse_part = initial::initial_partition(coarsest, nparts, &coarse_flat, options, &mut rng);

        // Refine at coarsest level using the configured algorithm
        refine_level(coarsest, &mut coarse_part, nparts, &coarse_flat, &ub, options, &SingleComm);

        let mut part = coarse_part;
        for i in (0..levels.len()).rev() {
            let level = &levels[i];
            let fine_graph = if i == 0 { graph } else { &levels[i - 1].graph };
            let fine_n = fine_graph.nvtxs;

            let mut fine_part = vec![0i32; fine_n];
            for v in 0..fine_n {
                fine_part[v] = part[level.cmap[v]];
            }

            let level_targets = fine_graph.compute_target_weights(nparts, tpwgts);
            let level_flat: Vec<f64> = level_targets.iter()
                .flat_map(|r| r.iter().copied()).collect();

            refine_level(fine_graph, &mut fine_part, nparts, &level_flat, &ub, options, &SingleComm);
            part = fine_part;
        }

        let objval = graph.edge_cut(&part);
        let result = PartitionResult { part, objval };
        if best_result.as_ref().map(|r| objval < r.objval).unwrap_or(true) {
            best_result = Some(result);
        }
    }

    // ── Parallel reduction: select globally best partition across all ranks ──
    if nprocs > 1 {
        let local_objval: Idx = best_result.as_ref().map(|r| r.objval).unwrap_or(Idx::MAX);
        let global_best_objval = comm.all_reduce_min_i32(local_objval);

        let local_part: Vec<Idx> = if best_result
            .as_ref()
            .map(|r| r.objval == global_best_objval)
            .unwrap_or(false)
        {
            best_result.as_ref().unwrap().part.clone()
        } else {
            vec![]
        };

        let gathered = comm.gather_i32_vec(0, &local_part);

        let mut winning_part = if comm.rank() == 0 {
            gathered
                .into_iter()
                .find(|v| !v.is_empty())
                .unwrap_or_else(|| {
                    best_result
                        .as_ref()
                        .map(|r| r.part.clone())
                        .unwrap_or_default()
                })
        } else {
            vec![]
        };

        comm.broadcast_i32_vec(0, &mut winning_part);

        return Ok(PartitionResult {
            objval: global_best_objval,
            part: winning_part,
        });
    }

    Ok(best_result.unwrap())
}

/// Apply the configured refinement algorithm for one multilevel level.
/// Greedy uses the distributed variant; FM falls back to the serial implementation.
fn refine_level(
    graph: &Graph,
    part: &mut Vec<Idx>,
    nparts: usize,
    tpwgts: &[f64],
    ubvec: &[f64],
    options: &Options,
    comm: &dyn Comm,
) {
    match options.rtype {
        RType::Greedy => {
            greedy_refine_with_comm(graph, part, nparts, tpwgts, ubvec, options.niter, comm);
        }
        RType::Fm => {
            refine::refine_kway(graph, part, nparts, tpwgts, ubvec, options);
        }
    }
}
