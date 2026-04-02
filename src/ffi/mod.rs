//! FFI modules.

#[cfg(feature = "c-api")]
pub mod c_api;

#[cfg(feature = "wasm")]
pub mod wasm_api;
