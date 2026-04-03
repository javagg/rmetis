//! MPI communication abstraction layer.
//!
//! Provides a [`Comm`] trait that wraps both the `jsmpi` (wasm32) and `mpi`
//! (native) backends behind a single interface. A [`SingleComm`] no-op
//! implementation is always available for single-process use.

use crate::types::Idx;

// ─── Public trait ─────────────────────────────────────────────────────────────

/// Abstract MPI communicator. All collective operations are blocking.
pub trait Comm: Send + Sync {
    /// This process's rank within the communicator.
    fn rank(&self) -> i32;

    /// Total number of processes.
    fn size(&self) -> i32;

    /// All-reduce with min operation; every rank gets the global minimum.
    fn all_reduce_min_i32(&self, local: Idx) -> Idx;

    /// All-reduce with sum over a flat `i32` slice; every rank gets the
    /// element-wise global sum. All ranks must pass equal-length slices.
    fn all_reduce_sum_i32_slice(&self, local: &[Idx], out: &mut [Idx]);

    /// Root broadcasts a `Vec<i32>` to all other ranks. Non-root ranks may
    /// pass an empty `Vec`; it is replaced by the broadcast data.
    fn broadcast_i32_vec(&self, root: i32, data: &mut Vec<Idx>);

    /// Each rank contributes one `Vec<i32>`. Root (rank 0) collects all
    /// contributions; non-root ranks get an empty outer `Vec`.
    fn gather_i32_vec(&self, root: i32, local: &[Idx]) -> Vec<Vec<Idx>>;

    /// All ranks contribute one `Vec<i32>`; every rank receives all contributions.
    fn all_gather_i32_vec(&self, local: &[Idx]) -> Vec<Vec<Idx>>;

    /// Global barrier.
    fn barrier(&self);
}

// ─── Single-process (no-op) ───────────────────────────────────────────────────

/// No-op communicator for single-process execution.
pub struct SingleComm;

impl Comm for SingleComm {
    fn rank(&self) -> i32 { 0 }
    fn size(&self) -> i32 { 1 }
    fn all_reduce_min_i32(&self, local: Idx) -> Idx { local }
    fn all_reduce_sum_i32_slice(&self, local: &[Idx], out: &mut [Idx]) {
        out.copy_from_slice(local);
    }
    fn broadcast_i32_vec(&self, _root: i32, _data: &mut Vec<Idx>) {}
    fn gather_i32_vec(&self, _root: i32, local: &[Idx]) -> Vec<Vec<Idx>> {
        vec![local.to_vec()]
    }
    fn all_gather_i32_vec(&self, local: &[Idx]) -> Vec<Vec<Idx>> {
        vec![local.to_vec()]
    }
    fn barrier(&self) {}
}

// ─── jsmpi backend (wasm32) ───────────────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
mod jsmpi_backend {
    use super::{Comm, Idx};
    use jsmpi::traits::{Communicator, Root};
    use jsmpi::SystemCommunicator;

    pub struct JsmpiComm {
        world: SystemCommunicator,
    }

    impl JsmpiComm {
        pub fn new(world: SystemCommunicator) -> Self {
            Self { world }
        }
    }

    impl Comm for JsmpiComm {
        fn rank(&self) -> i32 { self.world.rank() }
        fn size(&self) -> i32 { self.world.size() }

        fn all_reduce_min_i32(&self, local: Idx) -> Idx {
            // jsmpi supports only sum-reduce. Implement min via gather+broadcast.
            let gathered = self.gather_i32_vec(0, &[local]);
            let mut global_min = local;
            if self.rank() == 0 {
                global_min = gathered.into_iter()
                    .flat_map(|v| v)
                    .min()
                    .unwrap_or(local);
            }
            self.broadcast_i32_vec(0, &mut vec![global_min]);
            global_min
        }

        fn all_reduce_sum_i32_slice(&self, local: &[Idx], out: &mut [Idx]) {
            debug_assert_eq!(local.len(), out.len());
            let gathered = self.gather_i32_vec(0, local);
            let mut result: Vec<Idx> = vec![0; local.len()];
            if self.rank() == 0 {
                for row in &gathered {
                    for (i, &v) in row.iter().enumerate() {
                        result[i] += v;
                    }
                }
            }
            self.broadcast_i32_vec(0, &mut result);
            out.copy_from_slice(&result);
        }

        fn broadcast_i32_vec(&self, root_rank: i32, data: &mut Vec<Idx>) {
            let root = self.world.process_at_rank(root_rank);
            root.broadcast_into(data);
        }

        fn gather_i32_vec(&self, root_rank: i32, local: &[Idx]) -> Vec<Vec<Idx>> {
            let root = self.world.process_at_rank(root_rank);
            let local_vec: Vec<Idx> = local.to_vec();
            let mut out: Vec<Vec<Idx>> = Vec::new();
            root.gather_into_root(&local_vec, &mut out);
            out
        }

        fn all_gather_i32_vec(&self, local: &[Idx]) -> Vec<Vec<Idx>> {
            let mut result = self.gather_i32_vec(0, local);
            self.broadcast_i32_vec(0, &mut result.iter().flatten().cloned().collect());
            // Re-partition after broadcast
            let gathered = self.gather_i32_vec(0, local);
            let mut bcast = gathered;
            let root = self.world.process_at_rank(0);
            root.broadcast_into(&mut bcast);
            bcast
        }

        fn barrier(&self) {
            use jsmpi::traits::Communicator;
            self.world.barrier();
        }
    }
}

// ─── native mpi backend ──────────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
mod mpi_backend {
    use super::{Comm, Idx};
    use mpi::collective::CommunicatorCollectives;
    use mpi::collective::Root as MpiRoot;
    use mpi::collective::SystemOperation;
    use mpi::datatype::PartitionMut;
    use mpi::topology::Communicator as MpiCommunicator;
    use mpi::Count;

    pub struct MpiComm {
        universe: mpi::environment::Universe,
    }

    impl MpiComm {
        pub fn new(universe: mpi::environment::Universe) -> Self {
            Self { universe }
        }
    }

    impl Comm for MpiComm {
        fn rank(&self) -> i32 {
            self.universe.world().rank()
        }
        fn size(&self) -> i32 {
            self.universe.world().size()
        }

        fn all_reduce_min_i32(&self, local: Idx) -> Idx {
            let world = self.universe.world();
            let mut global = local;
            world.all_reduce_into(&local, &mut global, &SystemOperation::min());
            global
        }

        fn all_reduce_sum_i32_slice(&self, local: &[Idx], out: &mut [Idx]) {
            let world = self.universe.world();
            world.all_reduce_into(local, out, &SystemOperation::sum());
        }

        fn broadcast_i32_vec(&self, root_rank: i32, data: &mut Vec<Idx>) {
            let world = self.universe.world();
            let root = world.process_at_rank(root_rank);
            let mut len = data.len() as i32;
            root.broadcast_into(&mut len);
            if self.rank() != root_rank {
                data.resize(len as usize, 0);
            }
            root.broadcast_into(data.as_mut_slice());
        }

        fn gather_i32_vec(&self, root_rank: i32, local: &[Idx]) -> Vec<Vec<Idx>> {
            let world = self.universe.world();
            let root = world.process_at_rank(root_rank);
            let n = self.size() as usize;
            let local_len = local.len() as Count;

            if self.rank() == root_rank {
                // Root: gather lengths from all ranks
                let mut all_lens = vec![0 as Count; n];
                root.gather_into_root(&local_len, all_lens.as_mut_slice());

                let total: Count = all_lens.iter().sum();
                let mut buf = vec![0i32; total as usize];
                let displs: Vec<Count> = all_lens
                    .iter()
                    .scan(0 as Count, |acc, &l| {
                        let d = *acc;
                        *acc += l;
                        Some(d)
                    })
                    .collect();
                let mut partition = PartitionMut::new(buf.as_mut_slice(), all_lens.as_slice(), displs.as_slice());
                root.gather_varcount_into_root(local, &mut partition);

                let mut result = Vec::with_capacity(n);
                let mut offset = 0usize;
                for &l in &all_lens {
                    result.push(buf[offset..offset + l as usize].to_vec());
                    offset += l as usize;
                }
                result
            } else {
                root.gather_into(&local_len);
                root.gather_varcount_into(local);
                vec![]
            }
        }

        fn all_gather_i32_vec(&self, local: &[Idx]) -> Vec<Vec<Idx>> {
            let world = self.universe.world();
            let n = self.size() as usize;
            let local_len = local.len() as Count;

            // Allgather lengths
            let mut all_lens = vec![0 as Count; n];
            world.all_gather_into(&local_len, all_lens.as_mut_slice());

            let total: Count = all_lens.iter().sum();
            let mut buf = vec![0i32; total as usize];
            let displs: Vec<Count> = all_lens
                .iter()
                .scan(0 as Count, |acc, &l| {
                    let d = *acc;
                    *acc += l;
                    Some(d)
                })
                .collect();
            let mut partition = PartitionMut::new(buf.as_mut_slice(), all_lens.as_slice(), displs.as_slice());
            world.all_gather_varcount_into(local, &mut partition);

            let mut out = Vec::with_capacity(n);
            let mut offset = 0usize;
            for &l in &all_lens {
                out.push(buf[offset..offset + l as usize].to_vec());
                offset += l as usize;
            }
            out
        }

        fn barrier(&self) {
            self.universe.world().barrier();
        }
    }

    // Safety: Universe holds an MPI handle; MPI itself is thread-safe at the
    // MPI_THREAD_SERIALIZED level and above, which is what we need here.
    unsafe impl Send for MpiComm {}
    unsafe impl Sync for MpiComm {}
}

// ─── Public factory ──────────────────────────────────────────────────────────

/// Create a communicator wrapping the MPI world. Returns [`SingleComm`] when
/// MPI is absent or only one process is running.
pub fn world() -> Box<dyn Comm> {
    #[cfg(target_arch = "wasm32")]
    {
        match jsmpi::initialize() {
            Ok(universe) => {
                let w = universe.world();
                Box::new(jsmpi_backend::JsmpiComm::new(w))
            }
            Err(_) => Box::new(SingleComm),
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        use mpi::topology::Communicator as MpiCommunicator;
        match mpi::initialize() {
            Some(universe) => {
                if universe.world().size() > 1 {
                    Box::new(mpi_backend::MpiComm::new(universe))
                } else {
                    Box::new(SingleComm)
                }
            }
            None => Box::new(SingleComm),
        }
    }
}
