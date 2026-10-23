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

/// Allowed valences of the organic subset; the implicit hydrogen
/// count is the smallest allowed valence that absorbs the bond-order
/// sum, so hypervalent S/P chemistry (sulfate's S, phosphate's P)
/// needs no hydrogens where a single fixed valence would invent them.
pub fn allowed_valences(element: &str) -> &'static [f64] {
    match element {
        "C" => &[4.0],
        "N" => &[3.0],
        "O" => &[2.0],
        "S" => &[2.0, 4.0, 6.0],
        "P" => &[3.0, 5.0],
        "F" | "Cl" | "Br" | "I" | "H" => &[1.0],
        _ => &[],
    }
}

/// The valence-gap hydrogen count for one atom.
pub fn valence_gap(element: &str, bond_order_sum: f64) -> u32 {
    for allowed in allowed_valences(element) {
        if *allowed + 1e-9 >= bond_order_sum {
            return (allowed - bond_order_sum).round().max(0.0) as u32;
        }
    }
    0
}

/// The largest allowed valence.
pub fn valence(element: &str) -> f64 {
    allowed_valences(element).last().copied().unwrap_or(0.0)
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

    /// Number of incident bonds per atom (heavy-atom degree).
    pub fn degrees(&self) -> Vec<u32> {
        let mut degrees = vec![0; self.atoms.len()];
        for &(from, to, _) in &self.bonds {
            degrees[from] += 1;
            degrees[to] += 1;
        }
        degrees
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

    /// Implicit/explicit hydrogen count per atom: MDL H-count codes
    /// win, otherwise the valence gap (the oracle does the same for
    /// code-0 atoms — a bare `O` line derives H2O). Anionic
    /// heteroatoms take zero: their standardization crosses the
    /// `/p` layer.
    pub fn hydrogens(&self) -> Vec<u32> {
        let sums = self.bond_order_sums();
        self.atoms
            .iter()
            .zip(&sums)
            .map(|(atom, sum)| {
                if atom.h_count > 0 {
                    atom.h_count
                } else if atom.charge < 0 {
                    0
                } else {
                    valence_gap(&atom.element, *sum)
                }
            })
            .collect()
    }

    /// Element → count, with the hydrogens from [`Self::hydrogens`].
    pub fn composition(&self) -> Vec<(String, u32)> {
        use std::collections::BTreeMap;
        let mut map: BTreeMap<String, u32> = BTreeMap::new();
        for (atom, hydrogens) in self.atoms.iter().zip(self.hydrogens()) {
            *map.entry(atom.element.clone()).or_insert(0) += 1;
            if hydrogens > 0 {
                *map.entry("H".to_string()).or_insert(0) += hydrogens;
            }
        }
        map.into_iter().collect()
    }

    pub fn element(&self, index: usize) -> Result<&str> {
        self.atoms.get(index).map_or_else(
            || err(format!("atom index {index} out of range")),
            |a| Ok(a.element.as_str()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure { atoms, bonds }
    }

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            charge: 0,
        }
    }

    #[test]
    fn bare_oxygen_composes_water_the_way_the_oracle_does() {
        // A code-0 atom infers H from the valence gap: O with no
        // bonds -> 2 H (recorded identity/breadth water is H2O).
        let s = structure(vec![atom("O")], vec![]);
        assert_eq!(s.composition(), vec![("H".into(), 2), ("O".into(), 1)]);
    }

    #[test]
    fn mdl_hcount_codes_win_over_the_valence_gap() {
        // Ammonium: code 5 -> 4 explicit H (recorded as H3N only
        // after the proton standardization).
        let mut n = atom("N");
        n.h_count = 4;
        n.charge = 1;
        let s = structure(vec![n], vec![]);
        assert_eq!(s.composition(), vec![("H".into(), 4), ("N".into(), 1)]);
    }

    #[test]
    fn anionic_heteroatoms_take_no_implicit_hydrogen() {
        // Acetate: the anionic O carries no H; the mobile H enters via
        // the /p standardization, not here.
        let mut o_anion = atom("O");
        o_anion.charge = -1;
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), o_anion],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let composition = s.composition();
        assert_eq!(
            composition,
            vec![("C".into(), 2), ("H".into(), 3), ("O".into(), 2)],
            "acetate composition is C2H3O2 before standardization"
        );
    }

    #[test]
    fn aromatic_bonds_count_one_point_five() {
        // Benzene: each C sees 1.5 + 1.5 = 3.0 -> one valence gap -> 1 H.
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds = (0..6).map(|i| (i, (i + 1) % 6, 4)).collect();
        let s = structure(atoms, bonds);
        let composition = s.composition();
        assert_eq!(composition, vec![("C".into(), 6), ("H".into(), 6)]);
    }
}
