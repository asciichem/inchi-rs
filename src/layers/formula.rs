//! The `/formula` layer (M1): Hill order with the InChI twist
//! (carbon first when present, then hydrogen, then alphabetical;
//! counts omitted when 1), with the anion/cation hydrogen
//! standardization applied via the shared rules.

use crate::error::{err, Result};
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::LayerEmitter;

pub struct FormulaLayer;

impl LayerEmitter for FormulaLayer {
    fn name(&self) -> &'static str {
        "formula"
    }

    fn segment_kind(&self) -> Option<char> {
        None
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        let mut composition = structure.composition()?;
        let charge = structure.total_charge();
        if charge.abs() > 1 {
            return err(format!(
                "formula layer: |charge| = {charge} beyond the ±1 corpus gate"
            ));
        }
        let adjustment = standardization::hydrogen_adjustment(structure);
        if adjustment != 0 {
            add_hydrogen(&mut composition, adjustment);
        }
        Ok(Some(hill_string(&composition)))
    }
}

fn add_hydrogen(composition: &mut Vec<(String, u32)>, delta: i32) {
    for entry in composition.iter_mut() {
        if entry.0 == "H" {
            entry.1 = (entry.1 as i64 + delta as i64).max(0) as u32;
            return;
        }
    }
    if delta > 0 {
        composition.push(("H".to_string(), delta as u32));
    }
}

// Hill order with the InChI twist: carbon first (when present), then
// hydrogen, then the rest alphabetically; counts omitted when 1.
// Without carbon: pure alphabetical.
fn hill_string(composition: &[(String, u32)]) -> String {
    let has_carbon = composition.iter().any(|(e, _)| e == "C");
    let rank = |element: &str| match (has_carbon, element) {
        (true, "C") => 0,
        (true, "H") => 1,
        (true, _) => 2,
        (false, _) => 0,
    };
    let mut entries: Vec<&(String, u32)> = composition.iter().collect();
    entries.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)).then(a.0.cmp(&b.0)));
    entries
        .iter()
        .filter(|(_, count)| *count > 0)
        .map(|(element, count)| {
            if *count > 1 {
                format!("{element}{count}")
            } else {
                element.clone()
            }
        })
        .collect()
}
