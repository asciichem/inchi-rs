//! inchi-rs — the Rust InChI engine, gated byte-for-byte by the
//! oracle corpus in [inchi-tests](https://github.com/asciichem/inchi-tests).
//!
//! Architecture: an [`Engine`] façade over encapsulated concerns —
//! `molfile` reads the input edge, `graph` owns molecular-graph
//! semantics, and `layers` holds one emitter per InChI layer
//! (open for extension, closed for modification).

pub mod engine;
mod error;
pub mod graph;
pub mod layers;
mod molfile;

pub use engine::Engine;
pub use error::{Error, Result};
pub use graph::{AtomRecord, Structure};

/// Parses the V2000 molfile subset into a [`Structure`].
pub fn parse_molfile(input: &str) -> Result<Structure> {
    molfile::parse(input)
}
