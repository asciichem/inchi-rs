//! The ±1 proton standardization rules shared by the layers the
//! rules couple: the formula layer adjusts its hydrogens, the
//! proton layer emits the `/p` segment, the charge layer emits
//! `/q`. Decoded empirically from the oracle battery + corpus
//! (2026-10-25): a cation is proton-related iff it carries a
//! removable hydrogen or is a valence-deficient N/P center (an
//! all-heavy onium keeps its charge in `/q`); an anionic O/S/Se/Te/N
//! is proton-related only when it joins a mobile-H group (the parent
//! stores the H in the group — a lone alkoxide keeps `/q-1`);
//! opposite charges annihilate across a connected component (nitro,
//! N-oxide, betaine, diazo) and vanish from every layer.

use crate::canonical::tautomer;
use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

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

/// The decoded charge plan: how the structure's charges split
/// between the proton balance (`p`) and the fixed-charge identity
/// (`q`), with which side drove the proton shift.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChargePlan {
    pub q: i32,
    pub p: i32,
    pub p_cation: bool,
    pub p_anion: bool,
}

fn cation_is_proton_related(atom: &AtomRecord, degree: u32, bond_sum: f64) -> bool {
    atom.h_count > 0 || proton_deficient_pnicogen(atom, degree, bond_sum)
}

pub fn charge_plan(structure: &Structure) -> ChargePlan {
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    let adjusted = adjusted_hydrogens(structure);
    let mask = tautomer::endpoint_mask(structure, &adjusted);
    let parent_anion = first_charged_heteroatom(structure, false);
    let mut plan = ChargePlan {
        q: 0,
        p: 0,
        p_cation: false,
        p_anion: false,
    };
    for (index, atom) in structure.atoms.iter().enumerate() {
        if atom.charge > 0 {
            if cation_is_proton_related(atom, degrees[index], sums[index]) {
                plan.p += 1;
                plan.p_cation = true;
            } else {
                plan.q += 1;
            }
        } else if atom.charge < 0 {
            let hetero = matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N");
            let parent_path = parent_anion == Some(index);
            if hetero && (parent_path || mask[index]) {
                plan.p -= 1;
                plan.p_anion = true;
            } else {
                plan.q -= 1;
            }
        }
    }
    plan
}

/// The reference's internal-salt annihilation: an anionic O/S/Se/Te
/// cancels with a FIXED (proton-unrelated) cation anywhere in its
/// connected component — nitro N⁺–O⁻ and N-oxide (bonded), betaine
/// and acetylcarnitine's quaternary-N/carboxylate pairs (3-4 bonds).
/// A proton-related cation never annihilates: glycine's NH3⁺/COO⁻
/// pair stays, its charges netting in the proton balance. Runs before
/// any layer, so the annihilated anion registers no riding hydrogen
/// and no mobile group.
pub fn normalize_internal_salts(structure: &mut Structure) {
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); structure.atoms.len()];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    let mut consumed = vec![false; structure.atoms.len()];
    let anions: Vec<usize> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| {
            atom.charge < 0 && matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N")
        })
        .map(|(index, _)| index)
        .collect();
    for anion in anions {
        if consumed[anion] {
            continue;
        }
        // Component-wide walk from the anion: the first unconsumed
        // fixed cation annihilates with it.
        let mut visited = vec![false; structure.atoms.len()];
        visited[anion] = true;
        let mut queue = vec![anion];
        let mut found = None;
        while let Some(site) = queue.pop() {
            let is_cation = structure.atoms[site].charge > 0
                && !consumed[site]
                && !cation_is_proton_related(&structure.atoms[site], degrees[site], sums[site]);
            if is_cation {
                found = Some(site);
                break;
            }
            for &next in &neighbors[site] {
                if !visited[next] {
                    visited[next] = true;
                    queue.push(next);
                }
            }
        }
        if let Some(site) = found {
            // The reference never re-derives H after moving a charge:
            // the neutralized anion stays at its anionic reading
            // (zero implicit H), the cation loses its + to the pair.
            structure.atoms[anion].charge = 0;
            structure.atoms[anion].h_pinned = true;
            structure.atoms[site].charge = 0;
            consumed[anion] = true;
            consumed[site] = true;
        }
    }
}

fn proton_deficient_pnicogen(atom: &AtomRecord, degree: u32, bond_sum: f64) -> bool {
    matches!(atom.element.as_str(), "N" | "P") && degree <= 3 && bond_sum <= 3.0
}

/// The proton-layer value for this structure: `Some("p±n")` when the
/// standardization moved protons, `None` when the charge is part of
/// the identity (quaternary N / lone alkoxide / carbanion).
pub fn proton_layer(structure: &Structure) -> Option<String> {
    let p = charge_plan(structure).p;
    (p != 0).then(|| format!("p{p:+}"))
}

/// The charge-layer value: the fixed charges the normalization did
/// not move, in total.
pub fn charge_layer(structure: &Structure) -> Option<String> {
    let q = charge_plan(structure).q;
    (q != 0).then(|| format!("q{q:+}"))
}

/// The hydrogen adjustment the formula layer applies: an anion gains
/// its mobile H, a proton-carrying cation loses one.
pub fn hydrogen_adjustment(structure: &Structure) -> i32 {
    let plan = charge_plan(structure);
    match (plan.p_cation, plan.p_anion) {
        (true, false) => -1,
        (false, true) => 1,
        _ => 0,
    }
}

/// The hydrogen count per atom with the ±1 standardization applied at
/// its site. Anionic O/N split two ways: carried by carbon they take
/// the salt path — a riding hydrogen that survives only when the
/// mobile-H machinery finds the anion a group (the parent stores the
/// H in the group; a lone alkoxide keeps none) — while carried by
/// anything else (hydroxide, nitrate) they take the parent path: the
/// riding hydrogen stays FIXED, never mobile. A proton-related cation
/// loses one. The `/h` layer reads these; the formula layer's total
/// is their sum.
pub fn adjusted_hydrogens(structure: &Structure) -> Vec<u32> {
    let mut counts = structure.hydrogens();
    let parent_anion = first_charged_heteroatom(structure, false);
    let salt_anion = first_charged_heteroatom(structure, true);
    if let Some(index) = parent_anion.or(salt_anion) {
        counts[index] += 1;
    }
    if let Some(index) = first_cation_site(structure) {
        counts[index] = counts[index].saturating_sub(1);
    }
    if let Some(index) = salt_anion {
        let mask = tautomer::endpoint_mask(structure, &counts);
        if !mask[index] {
            counts[index] -= 1;
        }
    }
    counts
}

/// The first anionic O/N, optionally restricted to those carried by
/// carbon (the salt candidates; the reference's GetSaltChargeType
/// anchor condition).
fn first_charged_heteroatom(structure: &Structure, carbon_carried: bool) -> Option<usize> {
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); structure.atoms.len()];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    structure
        .atoms
        .iter()
        .enumerate()
        .position(|(index, atom)| {
            if atom.charge >= 0 || !(atom.element == "O" || atom.element == "N") {
                return false;
            }
            let carried = neighbors[index]
                .iter()
                .any(|&n| structure.atoms[n].element == "C");
            carried == carbon_carried
        })
}

fn first_cation_site(structure: &Structure) -> Option<usize> {
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    structure
        .atoms
        .iter()
        .enumerate()
        .position(|(index, atom)| {
            atom.charge > 0
                && (atom.h_count > 0
                    || proton_deficient_pnicogen(atom, degrees[index], sums[index]))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, charge: i32, h_count: u32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            charge,
            h_pinned: false,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            stereo_bonds: Default::default(),
            positions: Default::default(),
        }
    }

    #[test]
    fn protonated_amines_are_proton_related() {
        // Oracle battery: [NH4+] and C[NH+](C)C carry /p+1; the /h
        // parent loses the riding H.
        let s = structure(vec![atom("N", 1, 4)], vec![]);
        assert_eq!(proton_layer(&s), Some("p+1".to_string()));
        assert_eq!(charge_layer(&s), None);
        assert_eq!(adjusted_hydrogens(&s), vec![3]);
    }

    #[test]
    fn quaternary_nitrogen_is_fixed() {
        // C[N+](C)(C)C: /q+1, no /p, no hydrogen moved.
        let s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        );
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), Some("q+1".to_string()));
    }

    #[test]
    fn sp3_pnicogen_cations_without_hydrogen_are_proton_related() {
        // Quinuclidinium-class N+ (three single bonds, no H): the
        // charge is the proton deficit (fixture identity/breadth-233).
        let s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        );
        assert_eq!(proton_layer(&s), Some("p+1".to_string()));
        assert_eq!(charge_layer(&s), None);
    }

    #[test]
    fn lone_alkoxide_is_fixed_and_keeps_no_riding_hydrogen() {
        // [O-]C: /q-1; the formula stays CH3O — the riding H exists
        // only inside a mobile group.
        let s = structure(vec![atom("O", -1, 0), atom("C", 0, 3)], vec![(0, 1, 1)]);
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), Some("q-1".to_string()));
        assert_eq!(hydrogen_adjustment(&s), 0);
        assert_eq!(adjusted_hydrogens(&s), vec![0, 3]);
    }

    #[test]
    fn carboxylate_anion_joins_the_group_and_shifts_the_proton() {
        // CC(=O)[O-]: /p-1, formula gains the parent H, the anionic O
        // keeps its riding H (the mobile group holds it).
        let s = structure(
            vec![
                atom("C", 0, 3),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        assert_eq!(proton_layer(&s), Some("p-1".to_string()));
        assert_eq!(charge_layer(&s), None);
        assert_eq!(hydrogen_adjustment(&s), 1);
    }

    #[test]
    fn zwitterionic_glycine_nets_to_zero() {
        // [NH3+]CC(=O)[O-]: p = +1 (N) - 1 (carboxylate) = 0 — no
        // charge layers, identical to neutral glycine.
        let s = structure(
            vec![
                atom("N", 1, 3),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)],
        );
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), None);
    }

    #[test]
    fn betaine_annihilates() {
        // C[N+](C)(C)CC(=O)[O-]: the quaternary N and the
        // carboxylate O (three bonds apart) cancel — no layers.
        let mut s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![
                (0, 1, 1),
                (0, 2, 1),
                (0, 3, 1),
                (0, 4, 1),
                (4, 5, 1),
                (5, 6, 2),
                (5, 7, 1),
            ],
        );
        normalize_internal_salts(&mut s);
        assert!(s.atoms.iter().all(|a| a.charge == 0));
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), None);
    }

    #[test]
    fn protonated_cation_never_annihilates() {
        // Glycine's N is proton-related: the pair stays and nets in
        // the proton balance instead.
        let mut s = structure(
            vec![
                atom("N", 1, 3),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)],
        );
        normalize_internal_salts(&mut s);
        assert_eq!(s.atoms[0].charge, 1);
        assert_eq!(s.atoms[4].charge, -1);
    }
}
