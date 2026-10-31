//! The AuxInfo emitter: the auxiliary information the reference
//! writes beside every identifier — the canonical numbering of the
//! input atoms (`/N:`), the constitutional equivalence classes
//! (`/E:`), the normalized charges (`/CRV:`), and the reversibility
//! echo of the ORIGINAL submitted structure (`/rA`, `/rB`, `/rC`).
//! The canonical fields read the normalized structure the layers
//! speak; the echo fields are a pure function of the parse, explicit
//! hydrogens and all (the reference's `WriteOrigAtoms` /
//! `WriteOrigBonds` / `WriteOrigCoord`, ichiprt1.c).

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::hydrogen::canonical_ranks_for;
use crate::layers::standardization;

/// The AuxInfo line for a submitted structure (`original`, exactly as
/// parsed) and its normalized form (`working`, the one the layers
/// read).
pub fn aux_info(original: &Structure, working: &Structure) -> Result<String> {
    standardization::ensure_typable_charges(working)?;
    let counts = standardization::adjusted_hydrogens(working);
    let ranks = canonical_ranks_for(working);
    let groups = crate::canonical::tautomer::mobile_groups(working, &counts, &ranks);

    // The normalization-type digit: +1 a mobile-H representation
    // exists (groups with real mobile hydrogen — charge-only c-point
    // registrations are not t-groups and count nothing), +1 a group
    // carries BOTH mobile hydrogen and charge (the salt-merged
    // (Hn±) super-group — the dissociated-acid tautomerism), +4 a
    // hard (DE)PROTONATION moved a proton.
    let tautomeric = groups.iter().any(|group| group.count > 0);
    let salt_merged = groups
        .iter()
        .any(|group| group.count > 0 && group.charge != 0);
    let norm = i32::from(tautomeric)
        + i32::from(salt_merged)
        + 4 * i32::from(working.protons_hard_moved);

    let (adjacency, invariants) = crate::canonical::graph_of(working);
    // The aux speaks the numbering and symmetry of the representation
    // the reference fills AuxInfo from: for a tautomeric structure the
    // NON-TAUT one (mobile H removed from the endpoints, no t-group
    // rows), otherwise the structure as the layers read it.
    let ranks_aux;
    let symmetry;
    if tautomeric {
        let mask = crate::canonical::tautomer::endpoint_mask(working, &counts);
        let notaut = crate::canonical::tautomer::without_mobile_h(&counts, &mask);
        let (ranks, classes) =
            crate::canonical::canonical_ranks_and_symmetry(&adjacency, &invariants, &notaut, &[]);
        ranks_aux = ranks;
        symmetry = classes;
    } else {
        let groups_rows = crate::canonical::group_membership(working, &counts);
        let (ranks, classes) = crate::canonical::canonical_ranks_and_symmetry(
            &adjacency,
            &invariants,
            &counts,
            &groups_rows,
        );
        ranks_aux = ranks;
        symmetry = classes;
    }

    let mut out = format!("AuxInfo=1/{norm}/N:");
    // `/N:` — the original atom numbers in canonical order.
    let mut order: Vec<usize> = (0..working.atoms.len()).collect();
    order.sort_by_key(|&atom| ranks_aux[atom]);
    let numbering: Vec<String> = order.iter().map(|&atom| (atom + 1).to_string()).collect();
    out.push_str(&numbering.join(","));

    // `/E:` — the constitutional equivalence classes: members by
    // canonical number, non-singleton classes only, classes ordered
    // by their minimum.
    let mut classes: std::collections::BTreeMap<u32, Vec<u32>> = Default::default();
    for (atom, class) in symmetry.iter().enumerate() {
        classes.entry(*class).or_default().push(ranks_aux[atom]);
    }
    let mut equivalents: Vec<Vec<u32>> = classes
        .into_values()
        .filter(|members| members.len() > 1)
        .map(|mut members| {
            members.sort_unstable();
            members
        })
        .collect();
    equivalents.sort();
    if !equivalents.is_empty() {
        out.push_str("/E:");
        let text: Vec<String> = equivalents
            .iter()
            .map(|members| format!("({})", members.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",")))
            .collect();
        out.push_str(&text.join(""));
    }

    // `/gE:` — equivalence of the registered mobile groups: groups
    // whose members fill the same symmetry classes are one class,
    // labelled by the groups' 1-based positions in the identifier.
    if groups.len() > 1 {
        let atom_of_rank = |rank: u32| ranks_aux.iter().position(|&r| r == rank);
        let mut signatures: std::collections::BTreeMap<Vec<u32>, Vec<usize>> = Default::default();
        for (index, group) in groups.iter().enumerate() {
            let mut signature: Vec<u32> = group
                .ranks
                .iter()
                .filter_map(|&rank| atom_of_rank(rank))
                .map(|atom| symmetry[atom])
                .collect();
            signature.sort_unstable();
            signatures.entry(signature).or_default().push(index);
        }
        let mut equivalents: Vec<Vec<u32>> = signatures
            .into_values()
            .filter(|members| members.len() > 1)
            .map(|mut members| {
                members.sort_unstable();
                members.into_iter().map(|i| i as u32 + 1).collect()
            })
            .collect();
        equivalents.sort();
        if !equivalents.is_empty() {
            out.push_str("/gE:");
            let text: Vec<String> = equivalents
                .iter()
                .map(|members| {
                    format!(
                        "({})",
                        members.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",")
                    )
                })
                .collect();
            out.push_str(&text.join(""));
        }
    }

    // `/CRV:` — normalized charges and unusual valences, ascending
    // canonical number: `<canon>[±charge][.valence]`.
    let bond_sums = working.bond_order_sums();
    let degrees = working.degrees();
    let mut flagged: Vec<(u32, i32, i32)> = working
        .atoms
        .iter()
        .enumerate()
        .map(|(atom, a)| {
            let el = crate::canonical::protons::el_number(&a.element).unwrap_or(0);
            let unusual = crate::canonical::protons::unusual_el_valence(
                el,
                a.charge,
                bond_sums[atom].round() as i32,
                counts[atom],
                degrees[atom],
            );
            (ranks_aux[atom], a.charge, unusual)
        })
        .filter(|&(_, charge, unusual)| charge != 0 || unusual != 0)
        .collect();
    flagged.sort_unstable();
    if !flagged.is_empty() {
        out.push_str("/CRV:");
        let text: Vec<String> = flagged
            .iter()
            .map(|(rank, charge, unusual)| {
                let mut entry = rank.to_string();
                if *charge < 0 {
                    entry.push('-');
                } else if *charge > 0 {
                    entry.push('+');
                }
                if *charge != 0 {
                    entry.push_str(&charge.abs().to_string());
                }
                if *unusual != 0 {
                    entry.push('.');
                    entry.push_str(&unusual.to_string());
                }
                entry
            })
            .collect();
        out.push_str(&text.join(","));
    }

    out.push_str(&reversibility_block(original));
    Ok(out)
}

/// The `/rA` `/rB` `/rC` echo of the submitted structure, in the
/// input's own atom order, explicit hydrogens included.
fn reversibility_block(original: &Structure) -> String {
    let n = original.atoms.len();
    // The atom-block prefix letter records the chiral flag; the reader
    // does not carry it yet, so the structure reads non-chiral.
    let mut out = format!("/rA:{n}n");
    for atom in &original.atoms {
        out.push_str(&atom.element);
        match atom.charge {
            0 => {}
            negative if negative < 0 => {
                out.push('-');
                if negative < -1 {
                    out.push_str(&negative.abs().to_string());
                }
            }
            positive => {
                out.push('+');
                if positive > 1 {
                    out.push_str(&positive.to_string());
                }
            }
        }
        let dotted = atom.element.chars().count() == 1;
        if atom.isotope != 0 {
            if dotted {
                out.push('.');
            }
            out.push('i');
            out.push_str(&atom.isotope.to_string());
        }
    }

    // `/rB` — one `;`-terminated entry per atom 2..n: the bonds to
    // lower-numbered neighbors, ascending, each a type char (wedge
    // perspective applied) and the neighbor's 1-based number.
    out.push_str("/rB:");
    for j in 1..n {
        let mut entries: Vec<(usize, char)> = Vec::new();
        for (index, &(from, to, kind)) in original.bonds.iter().enumerate() {
            let (end, other) = if from == j {
                (j, to)
            } else if to == j {
                (j, from)
            } else {
                continue;
            };
            if other >= j {
                continue;
            }
            let stereo = original.stereo_bonds.get(&index).copied().unwrap_or(0);
            entries.push((other, bond_char(kind, stereo, end < other)));
        }
        entries.sort_unstable();
        for (neighbor, character) in entries {
            out.push(character);
            out.push_str(&(neighbor + 1).to_string());
        }
        out.push(';');
    }

    // `/rC` — one `;`-terminated coordinate per atom, trailing zeros
    // trimmed, an all-zero coordinate an empty entry.
    out.push_str("/rC:");
    for (x, y) in &original.positions {
        if *x == 0.0 && *y == 0.0 {
            out.push(';');
        } else {
            out.push_str(&format!("{},{},0;", clean_coord(*x), clean_coord(*y)));
        }
    }
    out
}

/// The bond-type character: `s`/`d`/`t`/`a`, with the wedge codes
/// rendered from the perspective of the endpoint currently listed —
/// up-at-this-end `p`, down-at-this-end `n`, either `v` (the
/// reference's WriteOrigBonds switch).
fn bond_char(kind: u32, stereo: u8, this_end_first: bool) -> char {
    match kind {
        2 => 'd',
        3 => 't',
        4 => 'a',
        _ => match stereo {
            1 => {
                if this_end_first {
                    'p'
                } else {
                    'n'
                }
            }
            6 => {
                if this_end_first {
                    'n'
                } else {
                    'p'
                }
            }
            4 => 'v',
            _ => 's',
        },
    }
}

/// The reference's `CleanOrigCoord`: fixed four-decimal coordinates
/// with trailing zeros (and the trailing point) trimmed; a leading
/// zero after the sign goes too — the oracle writes `-.7`.
fn clean_coord(value: f64) -> String {
    let mut text = format!("{value:.4}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    if let Some(rest) = text.strip_prefix("-0.") {
        return format!("-.{rest}");
    }
    if let Some(rest) = text.strip_prefix("0.") {
        return rest.to_string();
    }
    text
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
            h_pinned: false,
            isotope: 0,
        }
    }

    const WATER_MOL: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn water_reads_aux_info_with_empty_equivalence() {
        let original = crate::molfile::parse(WATER_MOL).expect("parses");
        let mut working = original.clone();
        working.fold_explicit_hydrogens();
        let aux = aux_info(&original, &working).expect("aux");
        assert_eq!(aux, "AuxInfo=1/0/N:1/rA:1nO/rB:;/rC:12,12,0;");
    }

    #[test]
    fn coord_trimming_matches_the_reference() {
        assert_eq!(clean_coord(12.0), "12");
        assert_eq!(clean_coord(1.4), "1.4");
        assert_eq!(clean_coord(-0.7), "-.7");
        assert_eq!(clean_coord(0.0), "0");
    }

    #[test]
    fn bond_chars_render_the_wedge_perspective() {
        assert_eq!(bond_char(1, 0, true), 's');
        assert_eq!(bond_char(2, 0, true), 'd');
        assert_eq!(bond_char(4, 0, true), 'a');
        assert_eq!(bond_char(1, 1, true), 'p');
        assert_eq!(bond_char(1, 1, false), 'n');
        assert_eq!(bond_char(1, 6, true), 'n');
        assert_eq!(bond_char(1, 6, false), 'p');
    }

    #[test]
    fn acetate_reads_its_mobile_group_and_charge() {
        // The recorded oracle line for the acetate battery fixture:
        // the mobile group turns the normalization digit to 1, the
        // anionic oxygen shows in /CRV: by its canonical number.
        let o_anion = AtomRecord {
            element: "O".to_string(),
            h_count: 0,
            charge: -1,
            h_pinned: false,
            isotope: 0,
        };
        let mut working = Structure {
            atoms: vec![atom("C"), atom("C"), atom("O"), o_anion],
            bonds: vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
            ..Default::default()
        };
        working.fold_explicit_hydrogens();
        let original = working.clone();
        let aux = aux_info(&original, &working).expect("aux");
        assert!(aux.starts_with("AuxInfo=1/1/N:"), "{aux}");
        assert!(aux.contains("/CRV:"), "{aux}");
        assert!(aux.contains("(3,4)"), "{aux}");
    }
}
