//! The molecular graph the layers derive from: atoms (element,
//! explicit H count, charge) plus bond orders, with the graph-level
//! queries layers need. Implicit-hydrogen and valence conventions
//! live here because they are graph semantics, not layer semantics.

use crate::error::{err, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct AtomRecord {
    pub element: String,
    /// Explicit H count: MDL code n means n−1 hydrogens (0 = infer
    /// from the valence gap).
    pub h_count: u32,
    pub charge: i32,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Structure {
    pub atoms: Vec<AtomRecord>,
    /// V2000 bond types: 1 single, 2 double, 3 triple, 4 aromatic
    /// (aromatic counts as 1.5 in valence sums).
    pub bonds: Vec<(usize, usize, u32)>,
}

/// Standard valences of the organic subset — the M1/M2 scope. S and
/// P entries exist but their oxoacid chemistry is mobile-H (M3)
/// territory.
pub fn valence(element: &str) -> f64 {
    match element {
        "C" => 4.0,
        "N" => 3.0,
        "O" => 2.0,
        "S" => 6.0,
        "P" => 5.0,
        "F" | "Cl" | "Br" | "I" | "H" => 1.0,
        _ => 0.0,
    }
}

pub fn bond_order(kind: u32) -> f64 {
    match kind {
        2 => 2.0,
        3 => 3.0,
        4 => 1.5,
        _ => 1.0,
    }
}

impl Structure {
    pub fn total_charge(&self) -> i32 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// Bond-order sums per atom index.
    pub fn bond_order_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.atoms.len()];
        for &(from, to, kind) in &self.bonds {
            let order = bond_order(kind);
            sums[from] += order;
            sums[to] += order;
        }
        sums
    }

    /// Element → count, with implicit hydrogens: MDL H-count codes
    /// win, otherwise the valence gap (the oracle does the same for
    /// code-0 atoms — a bare `O` line derives H2O). Anionic
    /// heteroatoms take zero: their standardization crosses the
    /// `/p` layer.
    pub fn composition(&self) -> Result<Vec<(String, u32)>> {
        use std::collections::BTreeMap;
        let sums = self.bond_order_sums();
        let mut map: BTreeMap<String, u32> = BTreeMap::new();
        for (atom, sum) in self.atoms.iter().zip(&sums) {
            *map.entry(atom.element.clone()).or_insert(0) += 1;
            let hydrogens = if atom.h_count > 0 {
                atom.h_count
            } else if atom.charge < 0 {
                0
            } else {
                (valence(&atom.element) - sum).max(0.0) as u32
            };
            if hydrogens > 0 {
                *map.entry("H".to_string()).or_insert(0) += hydrogens;
            }
        }
        Ok(map.into_iter().collect())
    }

    pub fn element(&self, index: usize) -> Result<&str> {
        self.atoms.get(index).map_or_else(
            || err(format!("atom index {index} out of range")),
            |a| Ok(a.element.as_str()),
        )
    }
}
