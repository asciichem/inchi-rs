//! Mobile-H (tautomer) detection — M3. A group is a set of O/N
//! endpoints linked through the documented 1,3 and 1,5 alternating
//! single/double shift paths: the donor terminal of every single bond
//! must be shift-capable (an attached H or a negative charge; a nitro
//! nitrogen never donates — the aci form keeps its H fixed), and
//! aromatic ring bonds match either parity, which is how the reference
//! delocalizes a ring before applying its templates. A group's
//! hydrogens are the ±1-standardized counts of its endpoints — the
//! acetate O⁻'s mobile H arrives exactly this way.

use crate::graph::Structure;
use std::collections::{BTreeSet, HashMap};

/// One mobile-H group: the hydrogen count and the canonical ranks of
/// its endpoints (ascending).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MobileGroup {
    pub count: u32,
    pub ranks: Vec<u32>,
}

fn is_endpoint(structure: &Structure, atom: usize) -> bool {
    matches!(structure.atoms[atom].element.as_str(), "O" | "N")
}

fn capable(structure: &Structure, hydrogens: &[u32], atom: usize) -> bool {
    if is_nitro_nitrogen(structure, atom) {
        // A nitro nitrogen's hydrogens stay fixed (the oracle's aci
        // form keeps them out of the mobile group).
        return false;
    }
    hydrogens[atom] > 0 || structure.atoms[atom].charge < 0
}

/// An N bonded to two or more O — the nitro situation.
fn is_nitro_nitrogen(structure: &Structure, atom: usize) -> bool {
    structure.atoms[atom].element == "N"
        && neighbors_of(structure, atom)
            .into_iter()
            .filter(|&(neighbor, _)| structure.atoms[neighbor].element == "O")
            .count()
            >= 2
}

fn neighbors_of(structure: &Structure, atom: usize) -> Vec<(usize, u32)> {
    structure
        .bonds
        .iter()
        .filter_map(|&(from, to, kind)| {
            if from == atom {
                Some((to, kind))
            } else if to == atom {
                Some((from, kind))
            } else {
                None
            }
        })
        .collect()
}

pub fn mobile_groups(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<MobileGroup> {
    let n = structure.atoms.len();
    let aromatic = crate::canonical::aromatic::aromatic_bonds(structure);
    let bond_index: HashMap<(usize, usize), usize> = structure
        .bonds
        .iter()
        .enumerate()
        .map(|(index, &(from, to, _))| ((from.min(to), from.max(to)), index))
        .collect();
    let double_like = |a: usize, b: usize| -> bool {
        bond_index
            .get(&(a.min(b), a.max(b)))
            .map(|&bond| structure.bonds[bond].2 == 2 || aromatic.contains(&bond))
            .unwrap_or(false)
    };

    let mut linked: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];

    // 1,3 pairs through a common neighbor: the single side must be
    // shift-capable; the anchor may itself be an endpoint (the nitro
    // anion's oxygens pair through their nitrogen).
    for anchor in 0..n {
        let ends: Vec<usize> = neighbors_of(structure, anchor)
            .into_iter()
            .filter(|&(atom, _)| is_endpoint(structure, atom))
            .map(|(atom, _)| atom)
            .collect();
        for &x in &ends {
            for &y in &ends {
                if y == x {
                    continue;
                }
                let xd = double_like(x, anchor);
                let yd = double_like(y, anchor);
                if (xd && (yd || capable(structure, hydrogens, y)))
                    || (yd && capable(structure, hydrogens, x))
                {
                    linked[x].insert(y);
                    linked[y].insert(x);
                }
            }
        }
    }

    // 1,5 walks (donor-s-d-s-d from a capable donor, either
    // direction); interior atoms may be endpoints.
    for start in 0..n {
        if !is_endpoint(structure, start) {
            continue;
        }
        if capable(structure, hydrogens, start) {
            walk_path(
                structure,
                &bond_index,
                &aromatic,
                hydrogens,
                start,
                start,
                1,
                FarRequirement::None,
                2,
                &mut linked,
            );
            // The 1,5 template terminates on a double bond — a bare
            // single-terminated acceptor is not one (phenol OH's stay
            // fixed).
            walk_path(
                structure,
                &bond_index,
                &aromatic,
                hydrogens,
                start,
                start,
                1,
                FarRequirement::Double,
                4,
                &mut linked,
            );
        }
        walk_path(
            structure,
            &bond_index,
            &aromatic,
            hydrogens,
            start,
            start,
            2,
            FarRequirement::Capable,
            2,
            &mut linked,
        );
        walk_path(
            structure,
            &bond_index,
            &aromatic,
            hydrogens,
            start,
            start,
            2,
            FarRequirement::Capable,
            4,
            &mut linked,
        );
    }

    let mut component = vec![usize::MAX; n];
    let mut groups = Vec::new();
    for start in 0..n {
        if !is_endpoint(structure, start) || component[start] != usize::MAX {
            continue;
        }
        let mut members = Vec::new();
        let mut queue = vec![start];
        component[start] = start;
        while let Some(atom) = queue.pop() {
            members.push(atom);
            for &next in &linked[atom] {
                if component[next] == usize::MAX {
                    component[next] = start;
                    queue.push(next);
                }
            }
        }
        let count: u32 = members.iter().map(|&atom| hydrogens[atom]).sum();
        if members.len() < 2 || count == 0 {
            continue;
        }
        let mut group_ranks: Vec<u32> = members.iter().map(|&atom| ranks[atom]).collect();
        group_ranks.sort_unstable();
        groups.push(MobileGroup {
            count,
            ranks: group_ranks,
        });
    }
    groups.sort_by_key(|group| (group.ranks.len(), group.ranks[0]));
    groups
}

/// What the far terminal of a shift walk must look like.
#[derive(Clone, Copy, PartialEq)]
enum FarRequirement {
    None,
    Capable,
    Double,
}

#[allow(clippy::too_many_arguments)]
fn walk_path(
    structure: &Structure,
    bond_index: &HashMap<(usize, usize), usize>,
    aromatic: &BTreeSet<usize>,
    hydrogens: &[u32],
    atom: usize,
    start: usize,
    parity: u32,
    far: FarRequirement,
    remaining: usize,
    linked: &mut Vec<BTreeSet<usize>>,
) {
    if remaining == 0 {
        return;
    }
    for (next, kind) in neighbors_of(structure, atom) {
        let key = (atom.min(next), atom.max(next));
        let aromatic_bond = bond_index
            .get(&key)
            .map(|&bond| aromatic.contains(&bond))
            .unwrap_or(false);
        if kind != parity && !aromatic_bond {
            continue;
        }
        if remaining == 1 {
            let far_key = (atom.min(next), atom.max(next));
            let far_double = bond_index
                .get(&far_key)
                .map(|&bond| structure.bonds[bond].2 == 2 || aromatic.contains(&bond))
                .unwrap_or(false);
            let accepted = match far {
                FarRequirement::None => true,
                FarRequirement::Capable => capable(structure, hydrogens, next),
                FarRequirement::Double => far_double,
            };
            if next == start || !is_endpoint(structure, next) || !accepted {
                continue;
            }
            linked[start].insert(next);
            linked[next].insert(start);
            continue;
        }
        walk_path(
            structure,
            bond_index,
            aromatic,
            hydrogens,
            next,
            start,
            if parity == 1 { 2 } else { 1 },
            far,
            remaining - 1,
            linked,
        );
    }
}

/// The atoms whose hydrogens live in mobile groups (excluded from the
/// fixed-H entries).
pub fn endpoint_atoms(ranks: &[u32], groups: &[MobileGroup]) -> BTreeSet<usize> {
    let mut by_rank = std::collections::BTreeMap::new();
    for (atom, &rank) in ranks.iter().enumerate() {
        by_rank.insert(rank, atom);
    }
    let mut endpoints = BTreeSet::new();
    for group in groups {
        for &rank in &group.ranks {
            if let Some(&atom) = by_rank.get(&rank) {
                endpoints.insert(atom);
            }
        }
    }
    endpoints
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            charge: 0,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure { atoms, bonds }
    }

    #[test]
    fn carboxyl_pairs_the_two_oxygens() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
    }

    #[test]
    fn alcohol_oxygen_is_not_an_endpoint() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O")],
            vec![(0, 1, 1), (1, 2, 1)],
        );
        assert!(mobile_groups(&s, &[3, 2, 1], &[1, 2, 3]).is_empty());
    }

    #[test]
    fn amide_nitrogen_joins_the_carbonyl_oxygen() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("N")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 2], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 2);
        assert_eq!(groups[0].ranks, vec![3, 4]);
        let endpoints = endpoint_atoms(&[1, 2, 3, 4], &groups);
        assert_eq!(endpoints, BTreeSet::from([2, 3]));
    }

    #[test]
    fn the_nitro_anion_groups_its_oxygens_but_not_its_nitrogen() {
        // CH3-N(=O)-O(-): the oracle records (H,3,4) over the two
        // oxygens; the nitrogen's hydrogens stay fixed.
        let s = structure(
            vec![atom("C"), atom("N"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[3, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
    }

    #[test]
    fn aromatic_bonds_carry_the_shift_across_a_ring() {
        // 2-hydroxypyridine skeleton: the OH oxygen and the ring N are
        // tautomeric through the aromatic ring.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms[0] = atom("N");
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 1), // ring N - C1
                (1, 2, 2),
                (2, 3, 1),
                (3, 4, 2),
                (4, 5, 1),
                (5, 0, 2), // ring closes with a double to N
                (1, 6, 1), // exocyclic OH on C1
            ],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 0, 0, 0, 1], &[1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![1, 7]);
    }
}
