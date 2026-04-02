//! Core types for rmetis: index types, options, errors, and result structures.

/// Index type — i32 to match METIS C ABI
pub type Idx = i32;
/// Real-valued weight type
pub type Real = f32;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// Coarsening scheme
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CType {
    /// Random matching
    Rm,
    /// Sorted heavy-edge matching
    Shem,
}

/// Initial partitioning scheme
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IpType {
    /// Graph-growing (BFS-based seeding)
    Grow,
    /// Random assignment
    Random,
    /// Recursive bisection via rmetis
    MetisRb,
}

/// Refinement scheme
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RType {
    /// Fiduccia-Mattheyses
    Fm,
    /// Greedy boundary refinement
    Greedy,
}

/// Objective function for k-way partitioning
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjType {
    /// Minimize edge cut
    Cut,
    /// Minimize communication volume
    Vol,
}

/// All tuning parameters — mirrors METIS options array semantics.
#[derive(Debug, Clone)]
pub struct Options {
    pub ctype: CType,
    pub iptype: IpType,
    pub rtype: RType,
    pub objtype: ObjType,
    /// Disable 2-hop matching in coarsening
    pub no2hop: bool,
    /// Number of independent partitioning attempts (best is kept)
    pub ncuts: usize,
    /// Number of refinement iterations per multilevel level
    pub niter: usize,
    /// Random seed (0 = entropy-seeded)
    pub seed: u64,
    /// Minimize maximum partition connectivity (k-way only)
    pub minconn: bool,
    /// Enforce contiguous partitions (k-way only)
    pub contig: bool,
    /// Imbalance tolerance × 1000; e.g. 30 → 3% slack. Range [1, 1000].
    pub ufactor: u32,
    /// 0 = 0-based numbering, 1 = 1-based (C ABI layer only)
    pub numbering: u8,
    /// Debug verbosity bitmask (0 = silent)
    pub dbglvl: u32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            ctype: CType::Shem,
            iptype: IpType::Grow,
            rtype: RType::Fm,
            objtype: ObjType::Cut,
            no2hop: false,
            ncuts: 1,
            niter: 10,
            seed: 0,
            minconn: false,
            contig: false,
            ufactor: 30,
            numbering: 0,
            dbglvl: 0,
        }
    }
}

impl Options {
    /// Options tuned for recursive bisection (matches METIS defaults for RB)
    pub fn for_recursive() -> Self {
        Options {
            iptype: IpType::MetisRb,
            rtype: RType::Fm,
            ufactor: 1,
            ..Default::default()
        }
    }

    /// Options tuned for k-way partitioning
    pub fn for_kway() -> Self {
        Options {
            iptype: IpType::Grow,
            rtype: RType::Fm,
            ufactor: 30,
            ..Default::default()
        }
    }
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------

/// Output of a graph partitioning operation
#[derive(Debug, Clone)]
pub struct PartitionResult {
    /// Partition assignment per vertex: `part[v]` in `[0, nparts)`
    pub part: Vec<Idx>,
    /// Objective value (edge cut or communication volume)
    pub objval: Idx,
}

/// Output of nested dissection ordering
#[derive(Debug, Clone)]
pub struct NDResult {
    /// Fill-reducing permutation: `perm[v]` = new position of original vertex v
    pub perm: Vec<Idx>,
    /// Inverse permutation: `iperm[i]` = original vertex at new position i
    pub iperm: Vec<Idx>,
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
pub enum MetisError {
    #[error("input error: {0}")]
    Input(String),

    #[error("out of memory")]
    OutOfMemory,

    #[error("graph is disconnected but contiguous partitions required")]
    Disconnected,

    #[error("cannot satisfy balance: min imbalance={0:.4}, required={1:.4}")]
    UnbalancedPartition(f32, f32),

    #[error("internal error: {0}")]
    Internal(String),
}

/// METIS-compatible C return codes
pub const METIS_OK: i32 = 1;
pub const METIS_ERROR_INPUT: i32 = -2;
pub const METIS_ERROR_MEMORY: i32 = -3;
pub const METIS_ERROR: i32 = -4;

/// Size of the METIS options array
pub const METIS_NOPTIONS: usize = 40;

// Option array indices (mirrors metis.h)
pub const METIS_OPTION_PTYPE: usize = 0;
pub const METIS_OPTION_OBJTYPE: usize = 1;
pub const METIS_OPTION_CTYPE: usize = 2;
pub const METIS_OPTION_IPTYPE: usize = 3;
pub const METIS_OPTION_RTYPE: usize = 4;
pub const METIS_OPTION_DBGLVL: usize = 5;
pub const METIS_OPTION_NITER: usize = 6;
pub const METIS_OPTION_NCUTS: usize = 7;
pub const METIS_OPTION_SEED: usize = 8;
pub const METIS_OPTION_NO2HOP: usize = 9;
pub const METIS_OPTION_MINCONN: usize = 10;
pub const METIS_OPTION_CONTIG: usize = 11;
pub const METIS_OPTION_UFACTOR: usize = 16;
pub const METIS_OPTION_NUMBERING: usize = 17;

impl Options {
    /// Parse from a raw METIS options array (C ABI compatibility)
    pub fn from_raw(arr: &[Idx; METIS_NOPTIONS]) -> Self {
        let mut opts = Options::default();
        if arr[METIS_OPTION_CTYPE] == 0 {
            opts.ctype = CType::Rm;
        }
        if arr[METIS_OPTION_IPTYPE] == 1 {
            opts.iptype = IpType::Random;
        } else if arr[METIS_OPTION_IPTYPE] == 4 {
            opts.iptype = IpType::MetisRb;
        }
        if arr[METIS_OPTION_RTYPE] == 1 {
            opts.rtype = RType::Greedy;
        }
        if arr[METIS_OPTION_OBJTYPE] == 1 {
            opts.objtype = ObjType::Vol;
        }
        if arr[METIS_OPTION_NO2HOP] != 0 {
            opts.no2hop = true;
        }
        if arr[METIS_OPTION_NCUTS] > 0 {
            opts.ncuts = arr[METIS_OPTION_NCUTS] as usize;
        }
        if arr[METIS_OPTION_NITER] > 0 {
            opts.niter = arr[METIS_OPTION_NITER] as usize;
        }
        opts.seed = arr[METIS_OPTION_SEED].max(0) as u64;
        opts.minconn = arr[METIS_OPTION_MINCONN] != 0;
        opts.contig = arr[METIS_OPTION_CONTIG] != 0;
        if arr[METIS_OPTION_UFACTOR] > 0 {
            opts.ufactor = arr[METIS_OPTION_UFACTOR] as u32;
        }
        opts.numbering = arr[METIS_OPTION_NUMBERING] as u8;
        opts.dbglvl = arr[METIS_OPTION_DBGLVL] as u32;
        opts
    }

    /// Serialize to raw METIS options array
    pub fn to_raw(&self) -> [Idx; METIS_NOPTIONS] {
        let mut arr = [-1i32; METIS_NOPTIONS];
        arr[METIS_OPTION_CTYPE] = match self.ctype { CType::Rm => 0, CType::Shem => 1 };
        arr[METIS_OPTION_IPTYPE] = match self.iptype {
            IpType::Grow => 0, IpType::Random => 1, IpType::MetisRb => 4,
        };
        arr[METIS_OPTION_RTYPE] = match self.rtype { RType::Fm => 0, RType::Greedy => 1 };
        arr[METIS_OPTION_OBJTYPE] = match self.objtype { ObjType::Cut => 0, ObjType::Vol => 1 };
        arr[METIS_OPTION_NO2HOP] = self.no2hop as Idx;
        arr[METIS_OPTION_NCUTS] = self.ncuts as Idx;
        arr[METIS_OPTION_NITER] = self.niter as Idx;
        arr[METIS_OPTION_SEED] = self.seed as Idx;
        arr[METIS_OPTION_MINCONN] = self.minconn as Idx;
        arr[METIS_OPTION_CONTIG] = self.contig as Idx;
        arr[METIS_OPTION_UFACTOR] = self.ufactor as Idx;
        arr[METIS_OPTION_NUMBERING] = self.numbering as Idx;
        arr[METIS_OPTION_DBGLVL] = self.dbglvl as Idx;
        arr
    }
}
