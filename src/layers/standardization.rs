//! The ±1 proton standardization rules shared by the layers the
//! rules couple: the formula layer adjusts its hydrogens, the
//! proton layer emits the `/p` segment. Decoded empirically from
//! the oracle corpus.

use crate::error::{err, Result};
use crate::graph::{valence, Structure};

/// The shared gate message: |charge| > 1 structures need the M3
/// mobile-H rules and are beyond the M1 layers' authority.
pub const BEYOND_M1: &str = "beyond the ±1 corpus gate";

/// Rejects structures whose |charge| exceeds the M1 gate — the ONE
/// place the policy lives; formula and proton layers both call it so
/// neither can silently half-decide a multi-charge structure.
pub fn ensure_within_m1(structure: &Structure) -> Result<()> {
    let charge = structure.total_charge();
    if charge.abs() > 1 {
        return err(format!(
            "standardization: |charge| = {} {BEYOND_M1}",
            charge.abs()
        ));
    }
    Ok(())
}

/// The proton-layer value for this structure: `Some("+1")` /
/// `Some("-1")` when the standardization moved a proton, `None` when
/// the charge is part of the identity (quaternary N / carbanion).
pub fn proton_layer(structure: &Structure) -> Option<String> {
    match structure.total_charge() {
        0 => None,
        -1 => {
            if hetero_anion(structure) {
                Some("-1".to_string())
            } else {
                None
            }
        }
        1 => {
            if cation_carries_hydrogen(structure) {
                Some("+1".to_string())
            } else {
                None
            }
        }
        _ => None,
    }
}

/// The hydrogen adjustment the formula layer applies: an anion gains
/// its mobile H, a proton-carrying cation loses one.
pub fn hydrogen_adjustment(structure: &Structure) -> i32 {
    match structure.total_charge() {
        -1 if hetero_anion(structure) => 1,
        1 if cation_carries_hydrogen(structure) => -1,
        _ => 0,
    }
}

fn hetero_anion(structure: &Structure) -> bool {
    structure
        .atoms
        .iter()
        .any(|a| a.charge < 0 && (a.element == "O" || a.element == "N"))
}

fn cation_carries_hydrogen(structure: &Structure) -> bool {
    let sums = structure.bond_order_sums();
    structure.atoms.iter().zip(&sums).any(|(atom, sum)| {
        atom.charge > 0 && (atom.h_count > 0 || (valence(&atom.element) - sum).max(0.0) > 0.0)
    })
}
