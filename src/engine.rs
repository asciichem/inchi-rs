//! The engine façade: parses a molfile and runs the registered
//! layers. Closed for modification, open for extension — new layers
//! join the registry in `layers::default_layers` without this file
//! changing.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::{default_layers, LayerEmitter};
use crate::molfile;

pub struct Engine {
    structure: Structure,
    layers: Vec<Box<dyn LayerEmitter>>,
}

impl Engine {
    /// Parses the V2000 molfile subset and equips the default layer
    /// registry.
    pub fn from_molfile(input: &str) -> Result<Self> {
        Ok(Engine {
            structure: molfile::parse(input)?,
            layers: default_layers(),
        })
    }

    /// The derived segments, tagged by layer kind — everything the
    /// registered layers emit. Milestones add layers; the gates match
    /// tagged segments against the oracle's recording by kind.
    pub fn segments(&self) -> Result<Vec<(Option<char>, String)>> {
        let mut out = Vec::new();
        for layer in &self.layers {
            if let Some(segment) = layer.emit(&self.structure)? {
                out.push((layer.segment_kind(), segment));
            }
        }
        Ok(out)
    }
}
