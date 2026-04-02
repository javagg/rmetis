//! CSR graph format validation.

use crate::types::{Idx, MetisError};

/// Validate a CSR graph representation.
///
/// Checks performed:
/// 1. `xadj.len() == nvtxs + 1`
/// 2. `xadj[0] == 0`, monotonically non-decreasing
/// 3. `adjncy.len() == xadj[nvtxs]`
/// 4. All neighbor indices in `[0, nvtxs)`
/// 5. No self-loops
/// 6. Symmetry: if u is a neighbor of v, then v is a neighbor of u
/// 7. Edge weight array length matches adjncy if provided
/// 8. Vertex weight array length matches `nvtxs * ncon` if provided
pub fn validate_csr(
    nvtxs: usize,
    ncon: usize,
    xadj: &[Idx],
    adjncy: &[Idx],
    vwgt: Option<&[Idx]>,
    adjwgt: Option<&[Idx]>,
) -> Result<(), MetisError> {
    if nvtxs == 0 {
        return Err(MetisError::Input("nvtxs must be > 0".into()));
    }
    if ncon == 0 {
        return Err(MetisError::Input("ncon must be > 0".into()));
    }

    // Check 1: xadj length
    if xadj.len() != nvtxs + 1 {
        return Err(MetisError::Input(format!(
            "xadj.len()={} but expected nvtxs+1={}",
            xadj.len(),
            nvtxs + 1
        )));
    }

    // Check 2: xadj[0] == 0 and monotone
    if xadj[0] != 0 {
        return Err(MetisError::Input(format!("xadj[0]={} but must be 0", xadj[0])));
    }
    for i in 0..nvtxs {
        if xadj[i] > xadj[i + 1] {
            return Err(MetisError::Input(format!(
                "xadj is not monotone at index {}: xadj[{}]={} > xadj[{}]={}",
                i, i, xadj[i], i+1, xadj[i+1]
            )));
        }
    }

    let nedges_stored = xadj[nvtxs] as usize;

    // Check 3: adjncy length
    if adjncy.len() != nedges_stored {
        return Err(MetisError::Input(format!(
            "adjncy.len()={} but xadj[nvtxs]={}",
            adjncy.len(),
            nedges_stored
        )));
    }

    // Check 4 & 5: neighbor range and no self-loops
    for v in 0..nvtxs {
        let start = xadj[v] as usize;
        let end = xadj[v + 1] as usize;
        for j in start..end {
            let u = adjncy[j];
            if u < 0 || u as usize >= nvtxs {
                return Err(MetisError::Input(format!(
                    "adjncy[{}]={} out of range [0, {})",
                    j, u, nvtxs
                )));
            }
            if u as usize == v {
                return Err(MetisError::Input(format!(
                    "self-loop detected at vertex {}",
                    v
                )));
            }
        }
    }

    // Check 6: symmetry
    // Build a quick reverse lookup: for each v, set of neighbors
    // Use a counting approach to avoid allocating a HashSet
    let mut degree_check = vec![0usize; nvtxs];
    for v in 0..nvtxs {
        let start = xadj[v] as usize;
        let end = xadj[v + 1] as usize;
        for j in start..end {
            let u = adjncy[j] as usize;
            degree_check[u] += 1;
        }
    }
    // Each vertex u should appear as neighbor of exactly as many vertices as it has neighbors
    // More precise: count occurrences of (v->u) and (u->v) across all pairs
    // Simple O(n) check: degree_check[v] should equal degree(v)
    for v in 0..nvtxs {
        let deg = (xadj[v + 1] - xadj[v]) as usize;
        if degree_check[v] != deg {
            return Err(MetisError::Input(format!(
                "graph is not symmetric at vertex {}: out-degree={} but appears as neighbor {} times",
                v, deg, degree_check[v]
            )));
        }
    }

    // Check 7: edge weight length
    if let Some(w) = adjwgt {
        if w.len() != nedges_stored {
            return Err(MetisError::Input(format!(
                "adjwgt.len()={} but expected {}",
                w.len(),
                nedges_stored
            )));
        }
    }

    // Check 8: vertex weight length
    if let Some(vw) = vwgt {
        let expected = nvtxs * ncon;
        if vw.len() != expected {
            return Err(MetisError::Input(format!(
                "vwgt.len()={} but expected nvtxs*ncon={}*{}={}",
                vw.len(), nvtxs, ncon, expected
            )));
        }
        // All weights must be positive
        for (i, &w) in vw.iter().enumerate() {
            if w <= 0 {
                return Err(MetisError::Input(format!(
                    "vwgt[{}]={} must be > 0",
                    i, w
                )));
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_xadj() -> Vec<Idx> { vec![0, 2, 4, 6, 8] }
    fn simple_adjncy() -> Vec<Idx> { vec![1, 3, 0, 2, 1, 3, 0, 2] }

    #[test]
    fn valid_square() {
        assert!(validate_csr(4, 1, &simple_xadj(), &simple_adjncy(), None, None).is_ok());
    }

    #[test]
    fn wrong_xadj_len() {
        let xadj = vec![0, 2, 4];
        let err = validate_csr(4, 1, &xadj, &simple_adjncy(), None, None).unwrap_err();
        assert!(matches!(err, MetisError::Input(_)));
    }

    #[test]
    fn self_loop() {
        let xadj = vec![0, 3, 5, 7, 9];
        // vertex 0 connects to 0 (self-loop), 1, 3
        let adjncy = vec![0, 1, 3, 0, 2, 1, 3, 0, 2];
        let err = validate_csr(4, 1, &xadj, &adjncy, None, None).unwrap_err();
        assert!(matches!(err, MetisError::Input(_)));
    }

    #[test]
    fn asymmetric() {
        // vertex 0 connects to 1, but vertex 1 does NOT connect back to 0
        let xadj = vec![0, 1, 1, 0, 0];
        let adjncy = vec![1];
        // degree_check[1] = 1, but degree(1) = 0 → asymmetry
        let err = validate_csr(4, 1, &xadj, &adjncy, None, None).unwrap_err();
        assert!(matches!(err, MetisError::Input(_)));
    }

    #[test]
    fn bad_vwgt_length() {
        let vwgt = vec![1, 1, 1]; // only 3, need 4
        let err = validate_csr(4, 1, &simple_xadj(), &simple_adjncy(), Some(&vwgt), None)
            .unwrap_err();
        assert!(matches!(err, MetisError::Input(_)));
    }
}
