//! C ABI compatibility layer.
//!
//! Provides METIS 5.1.x-compatible C function signatures.
//! This module is compiled only when the `c-api` feature is enabled.

use crate::graph::Graph;
use crate::partition::{kway, nd, recursive};
use crate::types::{
    Idx, MetisError, Options, Real,
    METIS_ERROR, METIS_ERROR_INPUT, METIS_ERROR_MEMORY, METIS_NOPTIONS, METIS_OK,
};

/// Fill an options array with default values.
#[no_mangle]
pub unsafe extern "C" fn RMETIS_SetDefaultOptions(options: *mut Idx) -> i32 {
    if options.is_null() {
        return METIS_ERROR_INPUT;
    }
    let opts = Options::default().to_raw();
    let slice = core::slice::from_raw_parts_mut(options, METIS_NOPTIONS);
    slice.copy_from_slice(&opts);
    METIS_OK
}

/// Recursive bisection graph partitioning.
#[no_mangle]
pub unsafe extern "C" fn RMETIS_PartGraphRecursive(
    nvtxs:   *const Idx,
    ncon:    *const Idx,
    xadj:    *const Idx,
    adjncy:  *const Idx,
    vwgt:    *const Idx,
    vsize:   *const Idx,
    adjwgt:  *const Idx,
    nparts:  *const Idx,
    tpwgts:  *const Real,
    ubvec:   *const Real,
    options: *const Idx,
    objval:  *mut Idx,
    part:    *mut Idx,
) -> i32 {
    let res = std::panic::catch_unwind(|| {
        c_partition_impl(
            nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt,
            nparts, tpwgts, ubvec, options, objval, part, false,
        )
    });
    match res {
        Ok(code) => code,
        Err(_) => METIS_ERROR,
    }
}

/// K-way graph partitioning.
#[no_mangle]
pub unsafe extern "C" fn RMETIS_PartGraphKway(
    nvtxs:   *const Idx,
    ncon:    *const Idx,
    xadj:    *const Idx,
    adjncy:  *const Idx,
    vwgt:    *const Idx,
    vsize:   *const Idx,
    adjwgt:  *const Idx,
    nparts:  *const Idx,
    tpwgts:  *const Real,
    ubvec:   *const Real,
    options: *const Idx,
    objval:  *mut Idx,
    part:    *mut Idx,
) -> i32 {
    let res = std::panic::catch_unwind(|| {
        c_partition_impl(
            nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt,
            nparts, tpwgts, ubvec, options, objval, part, true,
        )
    });
    match res {
        Ok(code) => code,
        Err(_) => METIS_ERROR,
    }
}

/// Nested dissection ordering.
#[no_mangle]
pub unsafe extern "C" fn RMETIS_NodeND(
    nvtxs:   *const Idx,
    xadj:    *const Idx,
    adjncy:  *const Idx,
    vwgt:    *const Idx,
    options: *const Idx,
    perm:    *mut Idx,
    iperm:   *mut Idx,
) -> i32 {
    let res = std::panic::catch_unwind(|| {
        if nvtxs.is_null() || xadj.is_null() || adjncy.is_null() || perm.is_null() || iperm.is_null() {
            return METIS_ERROR_INPUT;
        }
        let n = *nvtxs as usize;
        if n == 0 { return METIS_ERROR_INPUT; }

        let xadj_s = core::slice::from_raw_parts(xadj, n + 1);
        let ne = xadj_s[n] as usize;
        let adjncy_s = core::slice::from_raw_parts(adjncy, ne);
        let vwgt_opt = if vwgt.is_null() { None } else { Some(core::slice::from_raw_parts(vwgt, n)) };

        let opts = parse_options(options, METIS_NOPTIONS);

        let graph = match Graph::new(
            n, 1,
            xadj_s.to_vec(), adjncy_s.to_vec(),
            vwgt_opt.map(|s| s.to_vec()), None, None,
        ) {
            Ok(g) => g,
            Err(_) => return METIS_ERROR_INPUT,
        };

        match nd::node_nd(&graph, &opts) {
            Ok(result) => {
                let perm_s = core::slice::from_raw_parts_mut(perm, n);
                let iperm_s = core::slice::from_raw_parts_mut(iperm, n);
                perm_s.copy_from_slice(&result.perm);
                iperm_s.copy_from_slice(&result.iperm);
                METIS_OK
            }
            Err(MetisError::Input(_)) => METIS_ERROR_INPUT,
            Err(MetisError::OutOfMemory) => METIS_ERROR_MEMORY,
            Err(_) => METIS_ERROR,
        }
    });
    match res {
        Ok(code) => code,
        Err(_) => METIS_ERROR,
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

unsafe fn c_partition_impl(
    nvtxs:   *const Idx,
    ncon:    *const Idx,
    xadj:    *const Idx,
    adjncy:  *const Idx,
    vwgt:    *const Idx,
    vsize:   *const Idx,
    adjwgt:  *const Idx,
    nparts:  *const Idx,
    tpwgts:  *const Real,
    ubvec:   *const Real,
    options: *const Idx,
    objval:  *mut Idx,
    part:    *mut Idx,
    kway:    bool,
) -> i32 {
    // Mandatory pointer checks
    if nvtxs.is_null() || ncon.is_null() || xadj.is_null()
        || adjncy.is_null() || nparts.is_null() || part.is_null()
    {
        return METIS_ERROR_INPUT;
    }

    let n = *nvtxs as usize;
    let nc = *ncon as usize;
    let np = *nparts as usize;

    if n == 0 || nc == 0 || np < 2 { return METIS_ERROR_INPUT; }

    let xadj_s = core::slice::from_raw_parts(xadj, n + 1);
    let ne = xadj_s[n] as usize;
    let adjncy_s = core::slice::from_raw_parts(adjncy, ne);

    let vwgt_opt  = if vwgt.is_null()   { None } else { Some(core::slice::from_raw_parts(vwgt, n * nc).to_vec()) };
    let vsize_opt = if vsize.is_null()  { None } else { Some(core::slice::from_raw_parts(vsize, n).to_vec()) };
    let adjwgt_opt= if adjwgt.is_null() { None } else { Some(core::slice::from_raw_parts(adjwgt, ne).to_vec()) };
    let tpwgts_opt= if tpwgts.is_null() { None } else { Some(core::slice::from_raw_parts(tpwgts, np * nc)) };
    let ubvec_opt = if ubvec.is_null()  { None } else { Some(core::slice::from_raw_parts(ubvec, nc)) };

    let opts = parse_options(options, METIS_NOPTIONS);

    let graph = match Graph::new(
        n, nc,
        xadj_s.to_vec(), adjncy_s.to_vec(),
        vwgt_opt, adjwgt_opt, vsize_opt,
    ) {
        Ok(g) => g,
        Err(_) => return METIS_ERROR_INPUT,
    };

    let result = if kway {
        kway::partition_kway(&graph, np, tpwgts_opt, ubvec_opt, &opts)
    } else {
        recursive::partition_recursive(&graph, np, tpwgts_opt, ubvec_opt, &opts)
    };

    match result {
        Ok(r) => {
            let part_s = core::slice::from_raw_parts_mut(part, n);
            part_s.copy_from_slice(&r.part);
            if !objval.is_null() {
                *objval = r.objval;
            }
            METIS_OK
        }
        Err(MetisError::Input(_))           => METIS_ERROR_INPUT,
        Err(MetisError::OutOfMemory)         => METIS_ERROR_MEMORY,
        Err(_)                               => METIS_ERROR,
    }
}

unsafe fn parse_options(options: *const Idx, len: usize) -> Options {
    if options.is_null() {
        return Options::default();
    }
    let s = core::slice::from_raw_parts(options, len);
    // Check if array is initialized (first element != -1 signals custom options)
    if s[0] == -1 {
        return Options::default();
    }
    let arr: [Idx; crate::types::METIS_NOPTIONS] = s.try_into().unwrap_or([0; crate::types::METIS_NOPTIONS]);
    Options::from_raw(&arr)
}
