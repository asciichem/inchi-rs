//! Mobile-H (tautomer) detection — M3 slice 2. A group is a set of
//! O/N endpoints linked through common neighbors where at least one
//! endpoint-to-neighbor bond is a double bond (the 1,3-shift pattern:
//! O=C-OH, amide O=C-NH, guanidine N-C(=N)-N); linked endpoints merge
//! into one component (sulfuric-acid-like chains through S/P would
//! qualify but S/P valence is its own pending model). A group's
//! hydrogens are the ±1-standardized counts of its endpoints — the
//! acetate O⁻'s mobile H arrives exactly this way.

use crate::graph::Structure;
use std::collections::BTreeSet;

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
    hydrogens[atom] > 0 || structure.atoms[atom].charge < 0
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

#[allow(clippy::too_many_arguments)]
fn walk(
    structure: &Structure,
    hydrogens: &[u32],
    atom: usize,
    start: usize,
    parity: u32,
    far_must_be_capable: bool,
    remaining: usize,
    linked: &mut Vec<BTreeSet<usize>>,
) {
    if remaining == 0 {
        return;
    }
    for (next, kind) in neighbors_of(structure, atom) {
        if kind != parity {
            continue;
        }
        if remaining == 1 {
            if next == start
                || !is_endpoint(structure, next)
                || (far_must_be_capable && !capable(structure, hydrogens, next))
            {
                continue;
            }
            linked[start].insert(next);
            linked[next].insert(start);
            continue;
        }
        if is_endpoint(structure, next) {
            continue;
        }
        walk(
            structure,
            hydrogens,
            next,
            start,
            if parity == 1 { 2 } else { 1 },
            far_must_be_capable,
            remaining - 1,
            linked,
        );
    }
}

pub fn mobile_groups(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<MobileGroup> {
    let n = structure.atoms.len();
    // Linkage: X ~ Y when (a) some NON-endpoint atom Z bonds both and
    // at least one of the two bonds is a double — the 1,3 shift — or
    // (b) X and Y are directly double-bonded. Endpoints never anchor a
    // 1,3 link: an NH2 beside a tautomeric N stays fixed (breadth-24).
    let mut linked: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for start in 0..n {
        if !is_endpoint(structure, start) {
            continue;
        }
        // The reference's two documented shifts, walked from either
        // end: 1,3 (donor-s-d) and 1,5 (donor-s-d-s-d). The donor
        // terminal must be shift-capable; the far terminal needs
        // capability exactly when its own bond is single (the
        // double-start walks); intermediates are never endpoints.
        if capable(structure, hydrogens, start) {
            walk(structure, hydrogens, start, start, 1, false, 2, &mut linked);
            walk(structure, hydrogens, start, start, 1, false, 4, &mut linked);
        }
        walk(structure, hydrogens, start, start, 2, true, 2, &mut linked);
        walk(structure, hydrogens, start, start, 2, true, 4, &mut linked);
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
        // Ethanol: no double bond anywhere, no group.
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
    fn aromatic_bonds_do_not_anchor_groups() {
        // Benzene carbons are not O/N endpoints anyway; an aromatic O
        // (furan-like) must not group through a 1.5 bond.
        let atoms = (0..5).map(|_| atom("C")).collect();
        let mut s = structure(atoms, (0..5).map(|i| (i, (i + 1) % 5, 4)).collect());
        s.atoms[0] = atom("O");
        assert!(mobile_groups(&s, &[0, 0, 0, 0, 0], &[1, 2, 3, 4, 5]).is_empty());
    }
}
