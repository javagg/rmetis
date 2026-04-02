//! Multilevel k-way partitioning.

use crate::coarsen;
use crate::graph::Graph;
use crate::initial;
use crate::refine;
use crate::types::{MetisError, Options, PartitionResult, Real};
use rand::rngs::SmallRng;
use rand::SeedableRng;

/// Partition `graph` into `nparts` partitions using the multilevel k-way scheme.
pub fn partition_kway(
    graph: &Graph,
    nparts: usize,
    tpwgts: Option<&[Real]>,
    ubvec: Option<&[Real]>,
    options: &Options,
) -> Result<PartitionResult, MetisError> {
    if nparts < 2 {
        return Err(MetisError::Input(format!("nparts must be >= 2, got {}", nparts)));
    }
    if nparts > graph.nvtxs {
        return Err(MetisError::Input(format!(
            "nparts={} > nvtxs={}", nparts, graph.nvtxs
        )));
    }

    let ncon = graph.ncon;

    let targets = graph.compute_target_weights(nparts, tpwgts);
    let flat_targets: Vec<f64> = targets.iter().flat_map(|row| row.iter().copied()).collect();
    let ub: Vec<f64> = if let Some(ub) = ubvec {
        ub.iter().map(|&x| x as f64).collect()
    } else {
        let default_ub = 1.0 + options.ufactor as f64 / 1000.0;
        vec![default_ub; ncon]
    };

    let mut best_result: Option<PartitionResult> = None;

    for cut_attempt in 0..options.ncuts {
        let seed = if options.seed == 0 {
            let mut rng = SmallRng::from_entropy();
            use rand::Rng;
            rng.gen::<u64>()
        } else {
            options.seed.wrapping_add(cut_attempt as u64)
        };
        let mut rng = SmallRng::seed_from_u64(seed);

        // 1. Build coarsening hierarchy.
        // levels[0] = coarsened from original, levels.last() = coarsest.
        // Each level.cmap maps the *previous* (finer) graph's vertices to this level's vertices.
        // Previous graph for levels[0] is `graph`; for levels[i] (i>0) is levels[i-1].graph.
        let levels = coarsen::coarsen(graph, nparts, options, &mut rng);

        // 2. Initial partition on coarsest graph
        let coarsest = levels.last().map(|l| &l.graph).unwrap_or(graph);
        let coarse_targets = coarsest.compute_target_weights(nparts, tpwgts);
        let coarse_flat: Vec<f64> = coarse_targets.iter()
            .flat_map(|r| r.iter().copied()).collect();
        let mut coarse_part = initial::initial_partition(coarsest, nparts, &coarse_flat, options, &mut rng);

        // 3. Refine at coarsest level (use coarse-graph targets, not original targets)
        refine::refine_kway(coarsest, &mut coarse_part, nparts, &coarse_flat, &ub, options);

        // 4. Uncoarsen: project from coarsest back to original, refining at each level.
        // Iterate levels in reverse: coarsest→finest.
        // levels[i].cmap[v_fine] = v_coarse in levels[i].graph
        // We need the *finer* graph at each step. The finer graph for levels[i] is:
        //   - levels[i-1].graph if i > 0
        //   - the original `graph` if i == 0
        let mut part = coarse_part;
        for i in (0..levels.len()).rev() {
            let level = &levels[i];
            let fine_graph = if i == 0 { graph } else { &levels[i - 1].graph };
            let fine_n = fine_graph.nvtxs;

            // Project partition: fine vertex v → part[cmap[v]]
            let mut fine_part = vec![0i32; fine_n];
            for v in 0..fine_n {
                fine_part[v] = part[level.cmap[v]];
            }

            // Compute targets for this level's graph
            let level_targets = fine_graph.compute_target_weights(nparts, tpwgts);
            let level_flat: Vec<f64> = level_targets.iter()
                .flat_map(|r| r.iter().copied()).collect();

            // Refine on the finer graph
            refine::refine_kway(fine_graph, &mut fine_part, nparts, &level_flat, &ub, options);
            part = fine_part;
        }

        let objval = graph.edge_cut(&part);
        let result = PartitionResult { part, objval };
        if best_result.as_ref().map(|r| objval < r.objval).unwrap_or(true) {
            best_result = Some(result);
        }
    }

    Ok(best_result.unwrap())
}
