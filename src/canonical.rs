//! Canonical ranking, stage 1: invariant construction and iterative
//! refinement (one-dimensional Weisfeiler-Lehman). The refinement
//! yields the coarsest partition the invariants can distinguish; a
//! discrete partition is already a total order, while symmetric
//! graphs keep non-singleton classes for the stage-2
//! individualization search.

use crate::graph::{bond_order, Structure};

/// The initial invariant vector per atom. Strategies are traits so
/// the oracle corpus can arbitrate between invariant schemes without
/// touching the refinement.
pub trait InvariantScheme {
    fn initial(&self, structure: &Structure) -> Vec<Vec<u32>>;
}

/// InChI-style M2 invariants: element, charge, hydrogen count,
/// degree, and bond-order sum. The oracle ranks atoms by chemistry
/// first; these are the chemistry inputs available before any
/// connectivity refinement.
pub struct DefaultInvariants;

impl InvariantScheme for DefaultInvariants {
    fn initial(&self, structure: &Structure) -> Vec<Vec<u32>> {
        let hydrogens = structure.hydrogens();
        let degrees = structure.degrees();
        let sums = structure.bond_order_sums();
        (0..structure.atoms.len())
            .map(|index| {
                let atom = &structure.atoms[index];
                vec![
                    element_rank(&atom.element),
                    (atom.charge + 8).max(0) as u32,
                    hydrogens[index],
                    degrees[index],
                    half_units(sums[index]),
                ]
            })
            .collect()
    }
}

/// An atom partition: class id per atom. Ids are canonical
/// (assigned by sorting the signatures), so they do not depend on
/// atom order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Classes {
    ids: Vec<usize>,
    count: usize,
}

impl Classes {
    pub fn ids(&self) -> &[usize] {
        &self.ids
    }

    pub fn count(&self) -> usize {
        self.count
    }

    /// True when every atom has its own class — the partition is a
    /// total order and no individualization is needed.
    pub fn is_discrete(&self) -> bool {
        self.count == self.ids.len()
    }

    /// Atom indices in a class, ascending.
    pub fn members(&self, class: usize) -> impl Iterator<Item = usize> + '_ {
        self.ids
            .iter()
            .enumerate()
            .filter(move |(_, id)| **id == class)
            .map(|(index, _)| index)
    }
}

/// Refines the initial partition until it stops splitting.
pub fn refine(structure: &Structure, scheme: &dyn InvariantScheme) -> Classes {
    let mut classes = assign(scheme.initial(structure));
    loop {
        let next = assign(neighbor_signatures(structure, classes.ids()));
        if next.count() == classes.count() {
            return next;
        }
        classes = next;
    }
}

/// One refinement round: each atom is described by its class and the
/// sorted multiset of (bond order, neighbour class) pairs.
fn neighbor_signatures(structure: &Structure, ids: &[usize]) -> Vec<Vec<u32>> {
    let mut neighbors: Vec<Vec<(u32, u32)>> = vec![Vec::new(); structure.atoms.len()];
    for &(from, to, kind) in &structure.bonds {
        let order = half_units(bond_order(kind));
        neighbors[from].push((order, ids[to] as u32));
        neighbors[to].push((order, ids[from] as u32));
    }
    neighbors
        .into_iter()
        .enumerate()
        .map(|(index, mut list)| {
            list.sort_unstable();
            let mut signature = Vec::with_capacity(2 + list.len() * 2);
            signature.push(ids[index] as u32);
            signature.push(list.len() as u32);
            for (order, class) in list {
                signature.push(order);
                signature.push(class);
            }
            signature
        })
        .collect()
}

/// Maps signatures to class ids by their sorted order.
fn assign(signatures: Vec<Vec<u32>>) -> Classes {
    let mut distinct: Vec<&[u32]> = signatures.iter().map(Vec::as_slice).collect();
    distinct.sort_unstable();
    distinct.dedup();
    let ids = signatures
        .iter()
        .map(|signature| {
            distinct
                .binary_search(&signature.as_slice())
                .expect("signature was collected")
        })
        .collect();
    Classes {
        ids,
        count: distinct.len(),
    }
}

/// Bond-order sums are integral halves; keep the signature integral.
fn half_units(order: f64) -> u32 {
    (order * 2.0).round() as u32
}

/// Atomic-number order for the elements the corpus exercises; other
/// symbols fall back to a byte-packed, offset key. Class ids only
/// order signatures, so the fallback cannot change a partition.
fn element_rank(element: &str) -> u32 {
    const TABLE: [&str; 36] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr",
    ];
    if let Some(position) = TABLE.iter().position(|symbol| *symbol == element) {
        return position as u32 + 1;
    }
    let mut packed = 0u32;
    for byte in element.bytes().take(4) {
        packed = (packed << 8) | byte as u32;
    }
    1000 + packed
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(element: &str) -> crate::graph::AtomRecord {
        crate::graph::AtomRecord {
            element: element.to_string(),
            h_count: 0,
            charge: 0,
        }
    }

    fn structure(
        atoms: Vec<crate::graph::AtomRecord>,
        bonds: Vec<(usize, usize, u32)>,
    ) -> Structure {
        Structure { atoms, bonds }
    }

    #[test]
    fn benzene_is_one_symmetric_class() {
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds = (0..6).map(|i| (i, (i + 1) % 6, 4)).collect();
        let classes = refine(&structure(atoms, bonds), &DefaultInvariants);
        assert_eq!(classes.count(), 1);
        assert_eq!(classes.ids(), &[0; 6]);
        assert!(!classes.is_discrete());
    }

    #[test]
    fn ethanol_refines_to_a_total_order() {
        let classes = refine(
            &structure(
                vec![atom("C"), atom("C"), atom("O")],
                vec![(0, 1, 1), (1, 2, 1)],
            ),
            &DefaultInvariants,
        );
        assert_eq!(classes.count(), 3);
        assert!(classes.is_discrete());
    }

    #[test]
    fn isopropanol_keeps_the_methyls_together() {
        // C0 central, C1/C2 methyl, O3 hydroxyl.
        let classes = refine(
            &structure(
                vec![atom("C"), atom("C"), atom("C"), atom("O")],
                vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
            ),
            &DefaultInvariants,
        );
        assert_eq!(classes.count(), 3);
        assert_eq!(classes.ids()[1], classes.ids()[2]);
        assert_ne!(classes.ids()[0], classes.ids()[1]);
        assert_ne!(classes.ids()[0], classes.ids()[3]);
    }

    #[test]
    fn acetate_separates_its_four_heavy_atoms() {
        // CH3-C(=O)-O(-): methyl, carbonyl C, carbonyl O, anionic O.
        let mut anion = atom("O");
        anion.charge = -1;
        let classes = refine(
            &structure(
                vec![atom("C"), atom("C"), atom("O"), anion],
                vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
            ),
            &DefaultInvariants,
        );
        assert_eq!(classes.count(), 4);
        assert!(classes.is_discrete());
    }

    #[test]
    fn refinement_is_independent_of_atom_order() {
        // Isopropanol with the atom order rotated: old atom i becomes
        // atom (i + 1) % 4.
        let classes = refine(
            &structure(
                vec![atom("C"), atom("C"), atom("C"), atom("O")],
                vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
            ),
            &DefaultInvariants,
        );
        let rotated = |index: usize| (index + 1) % 4;
        let permuted = refine(
            &structure(
                vec![atom("O"), atom("C"), atom("C"), atom("C")],
                vec![
                    (rotated(0), rotated(1), 1),
                    (rotated(0), rotated(2), 1),
                    (rotated(0), rotated(3), 1),
                ],
            ),
            &DefaultInvariants,
        );
        for left in 0..4 {
            for right in 0..4 {
                assert_eq!(
                    classes.ids()[left] == classes.ids()[right],
                    permuted.ids()[rotated(left)] == permuted.ids()[rotated(right)],
                    "partition changed for atoms {left}/{right}"
                );
            }
        }
    }

    #[test]
    fn members_lists_the_class_atoms() {
        let classes = refine(
            &structure(
                vec![atom("C"), atom("C"), atom("C"), atom("O")],
                vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
            ),
            &DefaultInvariants,
        );
        let methyl = classes.ids()[1];
        assert_eq!(classes.members(methyl).collect::<Vec<_>>(), vec![1, 2]);
    }
}
