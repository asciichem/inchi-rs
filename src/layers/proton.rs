//! The `/p` layer (M1 scope): the ±1 proton standardization marker.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::LayerEmitter;

pub struct ProtonLayer;

impl LayerEmitter for ProtonLayer {
    fn name(&self) -> &'static str {
        "proton"
    }

    fn segment_kind(&self) -> Option<char> {
        Some('p')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        // The segment is the full InChI text: `p+1` / `p-1`.
        Ok(standardization::proton_layer(structure).map(|sign| format!("p{sign}")))
    }
}
