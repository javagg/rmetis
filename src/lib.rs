//! rmetis — Pure Rust graph partitioning library.
//!
//! Implements METIS 5.1.x-compatible graph partitioning algorithms:
//! - [`part_graph_kway`] — Multilevel k-way partitioning
//! - [`part_graph_recursive`] — Multilevel recursive bisection
//! - [`node_nd`] — Nested dissection for fill-reducing matrix ordering
//!
//! # Quick start
//!
//! ```rust
//! use rmetis::{Graph, Options, part_graph_kway};
//!
//! // Build a 4-vertex cycle graph (0-1-2-3-0)
//! let graph = Graph::new_unweighted(
//!     4,
//!     vec![0, 2, 4, 6, 8],
//!     vec![1, 3, 0, 2, 1, 3, 0, 2],
//! ).unwrap();
//!
//! let result = part_graph_kway(&graph, 2, None, None, &Options::default()).unwrap();
//! println!("partition: {:?}", result.part);
//! println!("edge cut:  {}", result.objval);
//! ```

pub mod coarsen;
pub mod ffi;
pub mod graph;
pub mod initial;
pub mod partition;
pub mod refine;
pub mod separator;
pub mod types;

// Re-export primary types
pub use graph::Graph;
pub use types::{
    CType, IpType, MetisError, NDResult, ObjType, Options, PartitionResult, RType,
    Idx, Real,
};

// Re-export primary API functions
pub use partition::kway::partition_kway as part_graph_kway;
pub use partition::recursive::partition_recursive as part_graph_recursive;
pub use partition::nd::node_nd;
