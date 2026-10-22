//! The ±1 proton standardization rules shared by the layers the
//! rules couple: the formula layer adjusts its hydrogens, the
//! proton layer emits the `/p` segment. Decoded empirically from
//! the oracle corpus.

use crate::graph::Structure;

/// The proton-layer value for this structure: `Some("+1")` /
/// `Some("-1")` when the standardization moved a proton, `None` when
/// the charge is part of the identity (quaternary N / carbanion).
pub fn proton_layer(structure: &Structure) -> Option<String> {
    match total_charge(structure) {
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
        _ => None, // beyond ±1 is M3 scope; callers reject earlier
    }
}

/// The hydrogen adjustment the formula layer applies: an anion gains
/// its mobile H, a proton-carrying cation loses one.
pub fn hydrogen_adjustment(structure: &Structure) -> i32 {
    match total_charge(structure) {
        -1 if hetero_anion(structure) => 1,
        1 if cation_carries_hydrogen(structure) => -1,
        _ => 0,
    }
}

fn total_charge(structure: &Structure) -> i32 {
    structure.atoms.iter().map(|a| a.charge).sum()
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
        atom.charge > 0
            && (atom.h_count > 0 || (crate::graph::valence(&atom.element) - sum).max(0.0) > 0.0)
    })
}
