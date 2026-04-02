//! Multilevel recursive bisection.

use crate::coarsen;
use crate::graph::Graph;
use crate::initial;
use crate::refine;
use crate::types::{Idx, MetisError, Options, PartitionResult, Real};
use rand::rngs::SmallRng;
use rand::SeedableRng;

/// Partition `graph` into `nparts` partitions using recursive bisection.
pub fn partition_recursive(
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

    let seed = if options.seed == 0 {
        let mut rng = SmallRng::from_entropy();
        use rand::Rng;
        rng.gen::<u64>()
    } else {
        options.seed
    };
    let mut rng = SmallRng::seed_from_u64(seed);

    // Build flat tpwgts: [nparts * ncon]
    let ncon = graph.ncon;
    let flat_tpwgts: Vec<Real> = if let Some(tp) = tpwgts {
        tp.to_vec()
    } else {
        vec![1.0 / nparts as Real; nparts * ncon]
    };

    let ub: Vec<f64> = if let Some(ub) = ubvec {
        ub.iter().map(|&x| x as f64).collect()
    } else {
        let default_ub = 1.0 + options.ufactor as f64 / 1000.0;
        vec![default_ub; ncon]
    };

    let mut part = vec![0i32; graph.nvtxs];
    rb_recurse(graph, nparts, 0, &flat_tpwgts, &ub, options, &mut rng, &mut part)?;

    let objval = graph.edge_cut(&part);
    Ok(PartitionResult { part, objval })
}

/// Recursive bisection: partition vertices (identified by `part[v] == current_id`)
/// into `nparts` groups, using partition IDs starting at `part_offset`.
fn rb_recurse(
    graph: &Graph,
    nparts: usize,
    part_offset: usize,
    tpwgts: &[Real],       // [nparts * ncon], targets for this sub-problem
    ubvec: &[f64],
    options: &Options,
    rng: &mut SmallRng,
    part: &mut Vec<Idx>,   // global partition array for the original graph
) -> Result<(), MetisError> {
    if nparts == 1 {
        // Assign all vertices to part_offset
        for v in 0..graph.nvtxs {
            part[v] = part_offset as Idx;
        }
        return Ok(());
    }

    // If too few vertices to split, assign greedily by weight
    if graph.nvtxs <= nparts {
        let n = graph.nvtxs;
        // Assign one vertex per partition, remaining go to last partition
        for v in 0..n {
            part[v] = (part_offset + v.min(nparts - 1)) as Idx;
        }
        return Ok(());
    }

    // Split into k1 and k2
    let k1 = nparts / 2;
    let k2 = nparts - k1;
    let ncon = graph.ncon;

    // Target weights for the bisection: left side vs right side
    let left_weight: Vec<Real> = (0..ncon)
        .map(|c| tpwgts[..k1 * ncon].iter().skip(c).step_by(ncon).sum())
        .collect();
    let right_weight: Vec<Real> = (0..ncon)
        .map(|c| tpwgts[k1 * ncon..].iter().skip(c).step_by(ncon).sum())
        .collect();
    let total: Real = left_weight[0] + right_weight[0];

    let bisect_tpwgts: Vec<Real> = left_weight.iter().chain(right_weight.iter())
        .map(|&w| w / total)
        .collect();

    // Multilevel bisection of current graph
    let bisect_part = multilevel_bisect(graph, &bisect_tpwgts, ubvec, options, rng);

    // Extract subgraphs and recurse
    let (g0, _o2n0, n2o0) = graph.extract_subgraph(&bisect_part, 0);
    let (g1, _o2n1, n2o1) = graph.extract_subgraph(&bisect_part, 1);

    // Normalized tpwgts for left sub-problem
    let s0: Real = left_weight[0];
    let left_sub_tp: Vec<Real> = tpwgts[..k1 * ncon].iter().map(|&w| w / s0).collect();
    let s1: Real = right_weight[0];
    let right_sub_tp: Vec<Real> = tpwgts[k1 * ncon..].iter().map(|&w| w / s1).collect();

    // Temporary sub-part arrays
    let mut sub_part0 = vec![0i32; g0.nvtxs];
    let mut sub_part1 = vec![0i32; g1.nvtxs];

    rb_recurse(&g0, k1, 0, &left_sub_tp, ubvec, options, rng, &mut sub_part0)?;
    rb_recurse(&g1, k2, 0, &right_sub_tp, ubvec, options, rng, &mut sub_part1)?;

    // Write results back to global part array
    for (new_v, &orig_v) in n2o0.iter().enumerate() {
        part[orig_v] = sub_part0[new_v] + part_offset as Idx;
    }
    for (new_v, &orig_v) in n2o1.iter().enumerate() {
        part[orig_v] = sub_part1[new_v] + (part_offset + k1) as Idx;
    }

    Ok(())
}

/// Single multilevel bisection (2-way partition).
fn multilevel_bisect(
    graph: &Graph,
    tpwgts: &[Real],  // [2 * ncon]
    ubvec: &[f64],
    options: &Options,
    rng: &mut SmallRng,
) -> Vec<Idx> {
    // Trivial bisection for very small graphs
    if graph.nvtxs <= 2 {
        let part: Vec<Idx> = (0..graph.nvtxs).map(|v| v as Idx).collect();
        return part;
    }

    let ncon = graph.ncon;
    let levels = coarsen::coarsen(graph, 2, options, rng);
    let coarsest = levels.last().map(|l| &l.graph).unwrap_or(graph);

    // Initial bisection
    let coarse_targets = coarsest.compute_target_weights(2, Some(tpwgts));
    let coarse_flat: Vec<f64> = coarse_targets.iter()
        .flat_map(|r| r.iter().copied()).collect();
    let mut part = initial::initial_partition(coarsest, 2, &coarse_flat, options, rng);

    let flat_targets: Vec<f64> = graph.compute_target_weights(2, Some(tpwgts))
        .iter().flat_map(|r| r.iter().copied()).collect();

    // Refine at coarsest level
    refine::refine_2way(coarsest, &mut part, 2, &coarse_flat, ubvec, options);

    // Uncoarsen + refine
    for i in (0..levels.len()).rev() {
        let level = &levels[i];
        let fine_graph = if i == 0 { graph } else { &levels[i - 1].graph };
        let fine_n = fine_graph.nvtxs;

        let mut fine_part = vec![0i32; fine_n];
        for v in 0..fine_n {
            fine_part[v] = part[level.cmap[v]];
        }

        let level_targets = fine_graph.compute_target_weights(2, Some(tpwgts));
        let level_flat: Vec<f64> = level_targets.iter()
            .flat_map(|r| r.iter().copied()).collect();
        refine::refine_2way(fine_graph, &mut fine_part, 2, &level_flat, ubvec, options);
        part = fine_part;
    }

    part
}
