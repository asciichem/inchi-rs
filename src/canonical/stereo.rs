//! Tetrahedral stereo (`/t`, `/m`, `/s`) — M4 slice 1, decoded from
//! oracle probes (the dev-time binary; input contract: a wedge (1) or
//! dash (6) bond-stereo code whose FIRST atom is the stereocenter —
//! exactly the gem writer's shape; the atom-parity field is ignored
//! for 2D input, and degenerate coordinates void wedges entirely).
//!
//! Parity model, fitted on the probe grid and the alanine fixtures:
//! the input ligand list is the heavy neighbors by atom-index order
//! with the implicit hydrogen FIRST; a center's parity is the wedge
//! direction (dash = 1) XOR the permutation parity of that list
//! reordered rank-ascending (H last). Emission normalizes the first
//! defined center's `/t` char to `-`, other defined centers show
//! `p_i ⊕ p_first` (`-` same, `+` differ), unwedged candidate centers
//! show `?`; `/m` is `p_first ⊕ (more than one center)`; `/s1` rides
//! `/m`. Two same-orbit centers with opposite parity are the meso
//! case — the oracle emits a bare `/t` with neither `/m` nor `/s`.

use crate::graph::Structure;

/// One center's derived state.
struct Center {
    atom: usize,
    /// Heavy-neighbor atom indices in input order (atom-index order).
    neighbors: Vec<usize>,
    /// Where the implicit H sits in the input list (always first when
    /// present — a 4-heavy-atom center has none).
    hydrogen: Option<usize>,
    parity: Option<u32>,
}

impl Center {
    /// The input ligand list as canonical ranks (None = the H).
    fn ligand_ranks(&self, ranks: &[u32]) -> Vec<Option<u32>> {
        let mut list: Vec<Option<u32>> = self.neighbors.iter().map(|&a| Some(ranks[a])).collect();
        if let Some(position) = self.hydrogen {
            list.insert(position, None);
        }
        list
    }
}

/// The derived segments in `/t`, `/m`, `/s` order; empty when the
/// structure carries no wedge parity at all.
pub fn segments(
    structure: &Structure,
    hydrogens: &[u32],
    ranks: &[u32],
    symmetry: &[u32],
) -> Vec<(char, String)> {
    let wedged = wedged_centers(structure);
    if wedged.is_empty() {
        return Vec::new();
    }
    let mut centers = candidates(structure, hydrogens, ranks);
    for &(atom, down) in &wedged {
        if let Some(center) = centers.iter_mut().find(|c| c.atom == atom) {
            let parity = permutation_parity(center.ligand_ranks(ranks));
            center.parity = Some(parity ^ down);
        }
    }
    let Some(first) = centers.iter().find(|c| c.parity.is_some()) else {
        return Vec::new();
    };
    let first_parity = first.parity.expect("checked");

    // Meso: same-orbit defined centers with opposite parity — the
    // relative case takes neither /m nor /s.
    let mut meso = false;
    for (i, a) in centers.iter().enumerate() {
        for b in centers.iter().skip(i + 1) {
            if let (Some(pa), Some(pb)) = (a.parity, b.parity) {
                if symmetry[a.atom] == symmetry[b.atom] && pa != pb {
                    meso = true;
                }
            }
        }
    }

    let mut t = String::from("t");
    for (index, center) in centers.iter().enumerate() {
        let mark = match center.parity {
            None => '?',
            Some(p) if index > 0 && p != first_parity => '+',
            Some(_) => '-',
        };
        if index > 0 {
            t.push(',');
        }
        t.push_str(&format!("{}{}", ranks[center.atom], mark));
    }
    let mut out = vec![('t', t)];
    if !meso {
        let flag = first_parity ^ u32::from(centers.len() > 1);
        out.push(('m', format!("m{flag}")));
        out.push(('s', "s1".to_string()));
    }
    out
}

/// Permutation parity (1 = odd) of the input list reordered
/// rank-ascending with the H (None) last.
fn permutation_parity(list: Vec<Option<u32>>) -> u32 {
    let mut remaining: Vec<(usize, Option<u32>)> = list.into_iter().enumerate().collect();
    remaining.sort_by_key(|&(position, rank)| (rank.is_none(), rank.unwrap_or(0), position));
    let sequence: Vec<usize> = remaining
        .into_iter()
        .map(|(position, _)| position)
        .collect();
    u32::try_from(inversions(&sequence) % 2).expect("parity fits u32")
}

fn inversions(sequence: &[usize]) -> usize {
    let mut count = 0;
    for i in 0..sequence.len() {
        for j in (i + 1)..sequence.len() {
            if sequence[i] > sequence[j] {
                count += 1;
            }
        }
    }
    count
}

/// Atoms the wedge codes point FROM (the code's first atom), with the
/// direction (1 = dash).
fn wedged_centers(structure: &Structure) -> Vec<(usize, u32)> {
    let mut out = Vec::new();
    for (&bond_index, &code) in &structure.stereo_bonds {
        let Some(&(from, _, _)) = structure.bonds.get(bond_index) else {
            continue;
        };
        let down = match code {
            6 => 1,
            1 => 0,
            _ => continue,
        };
        out.push((from, down));
    }
    out.sort();
    out
}

/// Every candidate center: three heavy neighbors plus one implicit H,
/// or four heavy neighbors, with all ligand ranks distinct.
fn candidates(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<Center> {
    let mut out = Vec::new();
    for (atom, &hydrogen_count) in hydrogens.iter().enumerate() {
        let neighbors: Vec<usize> = structure
            .bonds
            .iter()
            .filter_map(|&(from, to, _)| {
                if from == atom {
                    Some(to)
                } else if to == atom {
                    Some(from)
                } else {
                    None
                }
            })
            .collect();
        let hydrogen = if hydrogen_count == 1 {
            Some(0)
        } else if hydrogen_count == 0 && neighbors.len() == 4 {
            None
        } else {
            continue;
        };
        if hydrogen.is_some() && neighbors.len() != 3 {
            continue;
        }
        let mut ligands: Vec<u32> = neighbors.iter().map(|&a| ranks[a]).collect();
        ligands.sort_unstable();
        ligands.dedup();
        if ligands.len() != neighbors.len() {
            continue;
        }
        out.push(Center {
            atom,
            neighbors,
            hydrogen,
            parity: None,
        });
    }
    out.sort_by_key(|center| ranks[center.atom]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, h_count: u32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            charge: 0,
        }
    }

    fn butanol(wedge: u8) -> Structure {
        Structure {
            atoms: vec![
                atom("C", 3),
                atom("C", 1),
                atom("O", 2),
                atom("C", 2),
                atom("C", 3),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)],
            stereo_bonds: if wedge == 0 {
                Default::default()
            } else {
                [(1usize, wedge)].into_iter().collect()
            },
        }
    }

    #[test]
    fn the_probe_grid_reproduces_the_oracle_recordings() {
        // Canonical ranks from the recorded /c c1-3-4(2)5: the center
        // is rank 4, the wedged O ligand rank 2.
        let ranks = [1u32, 4, 2, 3, 5];
        let up = segments(&butanol(1), &[3, 1, 2, 2, 3], &ranks, &[1, 2, 3, 4, 5]);
        assert_eq!(
            up,
            vec![('t', "t4-".into()), ('m', "m1".into()), ('s', "s1".into())]
        );
        let down = segments(&butanol(6), &[3, 1, 2, 2, 3], &ranks, &[1, 2, 3, 4, 5]);
        assert_eq!(
            down,
            vec![('t', "t4-".into()), ('m', "m0".into()), ('s', "s1".into())]
        );
    }

    #[test]
    fn no_wedge_carries_no_stereo() {
        let ranks = [1u32, 4, 2, 3, 5];
        assert!(segments(&butanol(0), &[3, 1, 2, 2, 3], &ranks, &[1, 2, 3, 4, 5]).is_empty());
    }
}
