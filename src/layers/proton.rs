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

    fn order(&self) -> u8 {
        crate::layers::order::P
    }

    fn segment_kind(&self) -> Option<char> {
        Some('p')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_within_m1(structure)?;
        // The segment is the full InChI text: `p+1` / `p-1`.
        Ok(standardization::proton_layer(structure).map(|sign| format!("p{sign}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;
    use crate::layers::LayerEmitter;

    fn charged_atom(element: &str, charge: i32, h_count: u32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            charge,
        }
    }

    #[test]
    fn proton_layers_are_signed_segments() {
        let s = Structure {
            atoms: vec![charged_atom("N", 1, 4)],
            bonds: vec![],
        };
        assert_eq!(ProtonLayer.emit(&s).unwrap(), Some("p+1".to_string()));

        let s = Structure {
            atoms: vec![charged_atom("O", -1, 1)],
            bonds: vec![],
        };
        assert_eq!(ProtonLayer.emit(&s).unwrap(), Some("p-1".to_string()));
    }

    #[test]
    fn bond_charged_cations_emit_no_proton_layer() {
        // Quaternary N: four heavy bonds, no H — the + is the identity,
        // so no /p layer (fixture identity/tetramethylammonium).
        let s = Structure {
            atoms: vec![
                charged_atom("N", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        };
        assert_eq!(ProtonLayer.emit(&s).unwrap(), None);
    }
}
