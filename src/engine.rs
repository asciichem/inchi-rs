//! The engine façade: parses a molfile and runs the registered
//! layers. Closed for modification, open for extension — new layers
//! join the registry in `layers::default_layers` without this file
//! changing.

use crate::error::{err, Result};
use crate::graph::Structure;
use crate::layers::{default_layers, standardization, LayerEmitter, MANDATORY_KINDS};
use crate::molfile;

pub struct Engine {
    structure: Structure,
    layers: Vec<Box<dyn LayerEmitter>>,
}

impl Engine {
    /// Parses the V2000 molfile subset, normalizes internal salts,
    /// folds explicit hydrogen atoms onto their carriers, and equips
    /// the default layer registry.
    pub fn from_molfile(input: &str) -> Result<Self> {
        let mut structure = molfile::parse(input)?;
        standardization::normalize_internal_salts(&mut structure);
        structure.fold_explicit_hydrogens();
        Ok(Engine {
            structure,
            layers: default_layers(),
        })
    }

    /// The normalized structure — diagnostics for the triage tools.
    pub fn structure(&self) -> &Structure {
        &self.structure
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

    /// Assembles the full InChI string, segments sorted by declared
    /// layer order. Errors until the registered layers cover the
    /// mandatory kinds (formula + connection): a partial layer set
    /// must never present itself as an identifier.
    pub fn inchi_string(&self) -> Result<String> {
        let registered: Vec<Option<char>> = self.layers.iter().map(|l| l.segment_kind()).collect();
        for mandatory in &MANDATORY_KINDS {
            if !registered.contains(mandatory) {
                return err(format!(
                    "inchi_string: mandatory layer kind {:?} not registered yet — \
                     the engine's layer set is milestone-partial",
                    mandatory
                        .map(|c| c.to_string())
                        .unwrap_or_else(|| "formula".to_string())
                ));
            }
        }

        let mut sorted: Vec<&Box<dyn LayerEmitter>> = self.layers.iter().collect();
        sorted.sort_by_key(|layer| layer.order());
        let mut out = String::from("InChI=1S/");
        let mut appended = false;
        for layer in sorted {
            if let Some(segment) = layer.emit(&self.structure)? {
                if appended {
                    out.push('/');
                }
                out.push_str(&segment);
                appended = true;
            }
        }
        Ok(out)
    }

    /// The InChIKey over the assembled identifier — the truncated
    /// SHA-256 key, computed natively (see [`crate::key`]).
    pub fn inchikey(&self) -> Result<String> {
        crate::key::inchikey(&self.inchi_string()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Structure;
    use crate::layers::order;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    struct FixedSegment {
        order: u8,
        kind: Option<char>,
        text: &'static str,
    }

    impl LayerEmitter for FixedSegment {
        fn name(&self) -> &'static str {
            "fixed"
        }

        fn order(&self) -> u8 {
            self.order
        }

        fn segment_kind(&self) -> Option<char> {
            self.kind
        }

        fn emit(&self, _structure: &Structure) -> Result<Option<String>> {
            Ok(Some(self.text.to_string()))
        }
    }

    #[test]
    fn a_milestone_partial_layer_set_refuses_to_assemble_an_inchi() {
        // No connection layer registered: the set is milestone-partial
        // and must not present a prefix as an InChI.
        let engine = Engine {
            structure: crate::parse_molfile(WATER).expect("parses"),
            layers: vec![Box::new(FixedSegment {
                order: order::FORMULA,
                kind: None,
                text: "H2O",
            })],
        };
        let error = engine.inchi_string().unwrap_err();
        assert!(error.message().contains("mandatory layer kind"));
    }

    #[test]
    fn water_assembles_from_the_registered_layers() {
        let engine = Engine::from_molfile(WATER).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn assembly_follows_declared_layer_order_not_registry_order() {
        // Registry scrambled on purpose; the assembler must sort by
        // `order()` and join only the segments that emit.
        let engine = Engine {
            structure: Structure::default(),
            layers: vec![
                Box::new(FixedSegment {
                    order: order::P,
                    kind: Some('p'),
                    text: "p+1",
                }),
                Box::new(FixedSegment {
                    order: order::CONNECTION,
                    kind: Some('c'),
                    text: "c1-2",
                }),
                Box::new(FixedSegment {
                    order: order::FORMULA,
                    kind: None,
                    text: "H2O",
                }),
            ],
        };
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/c1-2/p+1");
    }
}
