//! The canonicalization search, mirroring the reference's two-stage
//! base flow.
//!
//! Stage 1 is hydrogen-blind: invariants, equitable refinement,
//! individualization, and the greatest connection table wins; every
//! table-equal terminal is an automorphic image of the first, and the
//! position-wise mappings merge into orbits. The orbit-minimal base
//! ranks become the first field of the stage-2 auxiliary ranks.
//!
//! Stage 2 rebuilds the initial partition from (base symmetry rank,
//! hydrogen count), refines it equitably, orders candidates within
//! cells by their auxiliary rank — the reference's `CellGetMinNode`
//! with `nAuxRank` — and keeps the greatest (table, hydrogens) terminal.
//! This stage is where hydrogen counts legitimately settle what pure
//! connectivity cannot.

use super::partition::Partition;
use super::table::ConnectionTable;

/// Orbit-minimal base ranks per atom: atoms the stage-1 search proved
/// automorphic share their smallest canonical rank — the symmetry
/// classes the stereo layer reads orbits from.
pub fn base_symmetry(adjacency: &[Vec<usize>], invariants: &[(u16, u16)]) -> Vec<u32> {
    let n = adjacency.len();
    let base = search_stage_one(adjacency, invariants);
    let mut symmetry: Vec<u32> = (0..n).map(|atom| base.partition.rank_of(atom)).collect();
    for atom in 0..n {
        let representative = base.orbits.representative(atom);
        symmetry[atom] = symmetry[representative].min(symmetry[atom]);
    }
    symmetry
}

/// Final canonical rank (1-based) per atom. The mobile groups (rank
/// free: count + member atoms) ride along as the terminal compare's
/// final rows — the reference's t-group CT section, which settles
/// the ties the connection table and the hydrogen layer leave equal.
pub fn canonical_ranks(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> Vec<u32> {
    let n = adjacency.len();
    let symmetry = base_symmetry(adjacency, invariants);
    // Stage-2 auxiliary partition: (symmetry rank, hydrogen count),
    // then equitable refinement — the reference's nRankAux.
    let aux_invariants: Vec<(u32, u32)> = (0..n)
        .map(|atom| (symmetry[atom], hydrogens[atom]))
        .collect();
    let mut aux = Partition::initial_pair(&aux_invariants);
    aux.refine(adjacency);
    let aux_rank: Vec<u32> = (0..n).map(|atom| aux.rank_of(atom)).collect();

    let informed = search_stage_two(adjacency, &aux, &aux_rank, hydrogens, groups);
    (0..n).map(|atom| informed.rank_of(atom)).collect()
}

/// The connection table over the final numbering.
pub fn canonical_table(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> ConnectionTable {
    let ranks = canonical_ranks(adjacency, invariants, hydrogens, groups);
    ConnectionTable::build(&ranks, adjacency)
}

/// Disjoint sets over atoms, merged by position-wise mappings of
/// automorphic terminal labelings.
struct Orbits {
    parent: Vec<usize>,
}

impl Orbits {
    fn new(n: usize) -> Orbits {
        Orbits {
            parent: (0..n).collect(),
        }
    }

    fn representative(&self, atom: usize) -> usize {
        let mut root = atom;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.representative(a), self.representative(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

struct StageOne {
    partition: Partition,
    orbits: Orbits,
}

fn search_stage_one(adjacency: &[Vec<usize>], invariants: &[(u16, u16)]) -> StageOne {
    let n = adjacency.len();
    let mut root = Partition::initial(invariants);
    root.refine(adjacency);
    let mut state = StageOneState {
        first: None,
        best: None,
        orbits: Orbits::new(n),
    };
    descend_hfree(&root, adjacency, &mut state);
    StageOne {
        partition: state.best.expect("a discrete partition is always reached"),
        orbits: state.orbits,
    }
}

struct StageOneState {
    first: Option<Partition>,
    best: Option<Partition>,
    orbits: Orbits,
}

fn descend_hfree(partition: &Partition, adjacency: &[Vec<usize>], state: &mut StageOneState) {
    match partition.first_nonsingleton_cell() {
        None => {
            let table = ConnectionTable::build(partition.ranks(), adjacency);
            if let Some(first) = &state.first {
                if table.entries() == ConnectionTable::build(first.ranks(), adjacency).entries() {
                    // Automorphic image of the first terminal: merge
                    // the position-wise mapping into the orbits.
                    for position in 0..adjacency.len() {
                        state
                            .orbits
                            .union(first.at_number()[position], partition.at_number()[position]);
                    }
                }
            } else {
                state.first = Some(partition.clone());
            }
            let better = match &state.best {
                None => true,
                Some(best_partition) => {
                    table.compare(&ConnectionTable::build(best_partition.ranks(), adjacency))
                        == std::cmp::Ordering::Greater
                }
            };
            if better {
                state.best = Some(partition.clone());
            }
        }
        Some((start, end)) => {
            for position in start..end {
                let v = partition.at_number()[position];
                let mut branch = partition.clone();
                branch.individualize(adjacency, v);
                descend_hfree(&branch, adjacency, state);
            }
        }
    }
}

fn search_stage_two(
    adjacency: &[Vec<usize>],
    aux: &Partition,
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> Partition {
    let n = adjacency.len();
    let mut best: Option<Partition> = None;
    let mut first: Option<Partition> = None;
    let mut orbits = Orbits::new(n);
    descend_hinformed(
        aux,
        adjacency,
        aux_rank,
        hydrogens,
        groups,
        &mut first,
        &mut best,
        &mut orbits,
    );
    best.expect("a discrete partition is always reached")
}

#[allow(clippy::too_many_arguments)]
fn descend_hinformed(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    first: &mut Option<Partition>,
    best: &mut Option<Partition>,
    orbits: &mut Orbits,
) {
    let cell = match partition.first_nonsingleton_cell() {
        Some(span) => span,
        None => {
            let ranks = partition.ranks();
            let table = ConnectionTable::build(ranks, adjacency);
            // Table-equal terminals are automorphic images: fold the
            // position-wise mapping into the orbit partition.
            let matches = |other: &Partition| {
                table.entries() == ConnectionTable::build(other.ranks(), adjacency).entries()
                    && (0..adjacency.len()).all(|position| {
                        hydrogens[other.at_number()[position]]
                            == hydrogens[partition.at_number()[position]]
                    })
            };
            if let Some(first_partition) = first {
                if matches(first_partition) {
                    for position in 0..adjacency.len() {
                        orbits.union(
                            first_partition.at_number()[position],
                            partition.at_number()[position],
                        );
                    }
                }
            } else {
                *first = Some(partition.clone());
            }
            let counts = |p: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|atom| hydrogens[p.at_number()[atom]])
                    .collect()
            };
            // The t-group CT rows: (count, ascending member ranks)
            // per group — the compare's final section.
            let rows = |p: &Partition| -> Vec<(u32, Vec<u32>)> {
                let ranks = p.ranks();
                groups
                    .iter()
                    .map(|(count, members)| {
                        let mut member_ranks: Vec<u32> =
                            members.iter().map(|&m| ranks[m]).collect();
                        member_ranks.sort_unstable();
                        (*count, member_ranks)
                    })
                    .collect()
            };
            let better = match &best {
                None => true,
                Some(best_partition) => {
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Less => false,
                        std::cmp::Ordering::Equal => {
                            match counts(partition).cmp(&counts(best_partition)) {
                                std::cmp::Ordering::Greater => true,
                                std::cmp::Ordering::Less => false,
                                std::cmp::Ordering::Equal => rows(partition) > rows(best_partition),
                            }
                        }
                    }
                }
            };
            if better {
                *best = Some(partition.clone());
            }
            return;
        }
    };
    // Candidates within the cell ordered by auxiliary rank, then atom
    // number — the reference's aux-informed CellGetMinNode — skipping
    // candidates an automorphism already relates to a tried one.
    let mut candidates: Vec<usize> = partition.at_number()[cell.0..cell.1].to_vec();
    candidates.sort_by_key(|&atom| (aux_rank[atom], atom));
    for &v in &candidates {
        if orbits.representative(v) != v {
            continue;
        }
        let mut branch = partition.clone();
        branch.individualize(adjacency, v);
        descend_hinformed(
            &branch, adjacency, aux_rank, hydrogens, groups, first, best, orbits,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::invariants::atom_invariants;

    fn run(elements: &[&str], bonds: &[(usize, usize)], hydrogens: &[u32]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        canonical_ranks(&adjacency, &invariants, hydrogens, &[])
    }

    fn table_of(elements: &[&str], bonds: &[(usize, usize)], hydrogens: &[u32]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        canonical_table(&adjacency, &invariants, hydrogens, &[])
            .entries()
            .to_vec()
    }

    #[test]
    fn acetaldehyde_numbers_methyl_then_carbonyl_chain() {
        // CH3-CH=O: the oracle records /c 1-2-3 with the methyl H3.
        assert_eq!(
            table_of(&["C", "C", "O"], &[(0, 1), (1, 2)], &[3, 1, 0]),
            vec![1, 2, 1, 3, 2]
        );
    }

    #[test]
    fn benzene_walks_two_then_four() {
        // The oracle records /c 1-2-4-6-5-3-1 around the ring.
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        assert_eq!(
            table_of(&["C"; 6], &bonds, &[1; 6]),
            vec![1, 2, 1, 3, 1, 4, 2, 5, 3, 6, 4, 5]
        );
    }

    #[test]
    fn isoprene_numbers_the_ch2s_before_the_methyl() {
        // The oracle records /c 1-4-5(2)3 with h4H,1-2H2,3H3: the
        // hydrogens order the hub's otherwise-equivalent neighbors.
        let ranks = run(
            &["C", "C", "C", "C", "C"],
            &[(0, 1), (1, 2), (3, 4), (1, 3)],
            &[3, 0, 2, 1, 2],
        );
        // The two CH2s are hydrogen-equivalent — either assignment
        // yields the recorded /c and /h.
        assert_eq!(ranks[0], 3);
        assert_eq!(ranks[1], 5);
        assert_eq!(ranks[3], 4);
        assert_eq!(sorted_pair(ranks[2], ranks[4]), (1, 2));
    }

    fn sorted_pair(a: u32, b: u32) -> (u32, u32) {
        if a <= b {
            (a, b)
        } else {
            (b, a)
        }
    }

    #[test]
    fn ranks_are_a_permutation() {
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        let mut ranks = run(&["C"; 6], &bonds, &[1; 6]);
        ranks.sort_unstable();
        assert_eq!(ranks, vec![1, 2, 3, 4, 5, 6]);
    }
}
