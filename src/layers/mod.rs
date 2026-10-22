//! Layer emitters — the engine's extension axis (OCP). Each InChI
//! layer is its own module implementing [`Layer`]; the engine runs
//! its registered layers without knowing any of their internals.
//! Adding a layer (connection, mobile-H, stereo) is a new module
//! plus one registration in [`default_layers`].

pub mod formula;
pub mod proton;
pub mod standardization;

use crate::error::Result;
use crate::graph::Structure;

pub trait LayerEmitter {
    /// The layer's name.
    fn name(&self) -> &'static str;

    /// The layer's segment identifier in the InChI string: `None`
    /// for the formula layer (the bare segment right after
    /// `InChI=1S/`), `Some('p')` for `/p`, `Some('c')` for `/c`, ...
    /// Gates match emitted segments to the oracle's recording BY
    /// THIS KIND, not by position — a milestone-partial engine emits
    /// a subset of the layers, and InChI orders them itself.
    fn segment_kind(&self) -> Option<char>;

    /// Emits the layer's segment, or `Ok(None)` when the structure
    /// carries nothing this layer records.
    fn emit(&self, structure: &Structure) -> Result<Option<String>>;
}

/// The engine's default layer set.
pub fn default_layers() -> Vec<Box<dyn LayerEmitter>> {
    vec![
        Box::new(formula::FormulaLayer),
        Box::new(proton::ProtonLayer),
    ]
}
