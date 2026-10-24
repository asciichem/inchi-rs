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
    /// The heavy ligands by ATOM NUMBER — the reference's neighbor[]
    /// order, which the half-parity is computed over.
    fn heavy_by_atom_number(&self) -> Vec<usize> {
        let mut sorted = self.neighbors.clone();
        sorted.sort_unstable();
        sorted
    }

    /// The permutation sequence the parity composes with: the heavy
    /// ligands' input positions reordered rank-ascending (the implicit
    /// H is rank-last and adds no transposition).
    fn permutation_sequence(&self, ranks: &[u32]) -> Vec<usize> {
        let mut remaining: Vec<(usize, u32)> = self
            .neighbors
            .iter()
            .enumerate()
            .map(|(position, &atom)| (position, ranks[atom]))
            .collect();
        remaining.sort_by_key(|&(_, rank)| rank);
        remaining
            .into_iter()
            .map(|(position, _)| position)
            .collect()
    }
}

/// The composition constant between the reference's half-bond parity
/// and our canonical parity, fitted on the alanine pair and the
/// 2-butanol grid.
const PARITY_CONSTANT: u32 = 0;

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
    for &(atom, bond_index, down) in &wedged {
        if let Some(center) = centers.iter_mut().find(|c| c.atom == atom) {
            // Heavy ligands by ATOM NUMBER (the reference's neighbor
            // order), lifted to 3D with the wedge out of the plane.
            let ligands = center.heavy_by_atom_number();
            let Some(half) = half_parity(structure, &ligands, bond_index, down) else {
                continue;
            };
            let perm = permutation_parity(center.permutation_sequence(ranks));
            center.parity = Some(half ^ perm ^ PARITY_CONSTANT);
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

/// Permutation parity (1 = odd) of the sequence that orders the heavy
/// ligands rank-ascending; the implicit H is rank-last and adds no
/// transposition.
fn permutation_parity(sequence: Vec<usize>) -> u32 {
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

/// The reference's half-stereo-bond parity (ichister.c
/// `half_stereo_bond_parity`): lift the three heavy ligand directions
/// to 3D — the wedge bond carries z = ± its drawn length, the
/// in-plane neighbors z = 0 — build a plane basis over the ligand
/// differences with the new z-axis on the world +z (wedge) side,
/// project, rotate ligand 0 onto +x, and read the counterclockwise
/// order of ligands 1 and 2 (CCW = odd = 1). `None` for degenerate
/// geometry (near-colinear bonds — the reference's ambiguity cases).
fn half_parity(
    structure: &Structure,
    ligands: &[usize],
    wedge_bond: usize,
    down: u32,
) -> Option<u32> {
    let position = |atom: usize| structure.positions.get(atom).copied();
    let (cx, cy) = position(structure.bonds[wedge_bond].0)?;
    let sign = if down == 1 { -1.0 } else { 1.0 };
    let mut coord = [[0.0f64; 3]; 3];
    for (slot, &atom) in ligands.iter().take(3).enumerate() {
        let (x, y) = position(atom)?;
        let (dx, dy) = (x - cx, y - cy);
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-9 {
            return None;
        }
        let z = if structure.bonds[wedge_bond].1 == atom {
            sign * len
        } else {
            0.0
        };
        coord[slot] = [dx / len, dy / len, z / len];
    }
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let len3 = |a: [f64; 3]| (dot(a, a)).sqrt();
    // Plane normal of the ligand-direction differences — kept as the
    // cross product gives it: the wedge's z sign lives in the normal's
    // x/y components, and forcing n.z > 0 would erase it.
    let d0 = sub(coord[0], coord[1]);
    let d1 = sub(coord[1], coord[2]);
    let mut normal = cross(d0, d1);
    let nlen = len3(normal);
    if nlen < 1e-9 {
        return None;
    }
    for value in &mut normal {
        *value /= nlen;
    }
    // Basis: ligand 0 projected onto the plane, then its right-handed
    // in-plane complement.
    let height = dot(coord[0], normal);
    let mut x_axis = sub(
        coord[0],
        [height * normal[0], height * normal[1], height * normal[2]],
    );
    let xlen = len3(x_axis);
    if xlen < 1e-9 {
        return None;
    }
    for value in &mut x_axis {
        *value /= xlen;
    }
    let y_axis = cross(normal, x_axis);
    // Ligand directions in the new basis, rotated so ligand 0 is +x.
    let project = |c: [f64; 3]| [dot(c, x_axis), dot(c, y_axis)];
    let p0 = project(coord[0]);
    let angle = |p: [f64; 2]| -> f64 {
        let c = p0[0] * p[0] + p0[1] * p[1];
        let s = -p0[1] * p[0] + p0[0] * p[1];
        let a = (s * p[0] + c * p[1]).atan2(c * p[0] + s * p[1]);
        if a < 0.0 {
            a + std::f64::consts::TAU
        } else {
            a
        }
    };
    let a1 = angle(project(coord[1]));
    let a2 = angle(project(coord[2]));
    let degenerate = (a1 - a2).abs().sin().abs() < 1e-9 && {
        let (min, max) = (a1.min(a2), a1.max(a2));
        min > std::f64::consts::PI - 1e-9 || max < std::f64::consts::PI + 1e-9
    };
    if degenerate {
        return None;
    }
    Some(u32::from(a1 < a2))
}

/// Atoms the wedge codes point FROM (the code's first atom), with the
/// bond index and the direction (1 = dash).
fn wedged_centers(structure: &Structure) -> Vec<(usize, usize, u32)> {
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
        out.push((from, bond_index, down));
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
            positions: vec![
                (12.0, 22.0),
                (72.0, 67.0),
                (132.0, 22.0),
                (72.0, 12.0),
                (132.0, 62.0),
            ],
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
