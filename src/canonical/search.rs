//! The individualization-refinement search: recursively split the
//! first non-singleton cell and keep the terminal partition with the
//! greatest connection table. The reference's automorphism pruning
//! only skips equal-table subtrees, so the winner is identical.
//!
//! Like the reference's `GetCanonRanking2`, the search runs twice:
//! the first pass fixes the atom order, the graph is renumbered into
//! it, and the second pass re-canonicalizes. The second pass is what
//! settles symmetric molecules — its exploration order differs from
//! the first pass's, and the corpus arbitrates the combination.

use super::partition::Partition;
use super::table::ConnectionTable;

/// Final canonical rank (1-based) per atom. `hydrogens` breaks the
/// ties connectivity cannot: among table-equal terminals the one
/// whose hydrogen counts sort lower wins (the corpus reading of the
/// reference's layered table comparison).
pub fn canonical_ranks(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
) -> Vec<u32> {
    let n = adjacency.len();
    let first = search_once(adjacency, invariants, hydrogens);
    let mut renumbered = vec![Vec::new(); n];
    for atom in 0..n {
        renumbered[first.rank_of(atom) as usize - 1] = adjacency[atom]
            .iter()
            .map(|&b| first.rank_of(b) as usize - 1)
            .collect();
    }
    let reordered: Vec<(u16, u16)> = first
        .at_number()
        .iter()
        .map(|&atom| invariants[atom])
        .collect();
    let reordered_hydrogens: Vec<u32> = first
        .at_number()
        .iter()
        .map(|&atom| hydrogens[atom])
        .collect();
    let second = search_once(&renumbered, &reordered, &reordered_hydrogens);
    (0..n)
        .map(|atom| second.rank_of(first.rank_of(atom) as usize - 1))
        .collect()
}

/// The connection table over the final numbering.
pub fn canonical_table(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
) -> ConnectionTable {
    let ranks = canonical_ranks(adjacency, invariants, hydrogens);
    ConnectionTable::build(&ranks, adjacency)
}

fn search_once(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
) -> Partition {
    let mut root = Partition::initial(invariants);
    root.refine(adjacency);
    let mut best: Option<Partition> = None;
    descend(&root, adjacency, hydrogens, &mut best);
    best.expect("a discrete partition is always reached")
}

fn descend(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    hydrogens: &[u32],
    best: &mut Option<Partition>,
) {
    match partition.first_nonsingleton_cell() {
        None => {
            let ranks = partition.ranks();
            let counts = |partition: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|atom| hydrogens[partition.at_number()[atom]])
                    .collect()
            };
            let better = match &best {
                None => true,
                Some(best_partition) => {
                    let table = ConnectionTable::build(ranks, adjacency);
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Equal => counts(partition) < counts(best_partition),
                        _ => false,
                    }
                }
            };
            if better {
                *best = Some(partition.clone());
            }
        }
        Some((start, end)) => {
            for position in start..end {
                let v = partition.at_number()[position];
                let mut branch = partition.clone();
                branch.individualize(adjacency, v);
                descend(&branch, adjacency, hydrogens, best);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::invariants::atom_invariants;

    fn ranks_of(elements: &[&str], bonds: &[(usize, usize)]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        let hydrogens = vec![0u32; elements.len()];
        canonical_ranks(&adjacency, &invariants, &hydrogens)
    }

    fn table_of(elements: &[&str], bonds: &[(usize, usize)]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        let hydrogens = vec![0u32; elements.len()];
        canonical_table(&adjacency, &invariants, &hydrogens)
            .entries()
            .to_vec()
    }

    #[test]
    fn acetaldehyde_numbers_methyl_then_carbonyl_chain() {
        // CH3-CH=O: the oracle records /c 1-2-3.
        assert_eq!(
            table_of(&["C", "C", "O"], &[(0, 1), (1, 2)]),
            vec![1, 2, 1, 3, 2]
        );
    }

    #[test]
    fn benzene_walks_two_then_four() {
        // The oracle records /c 1-2-4-6-5-3-1 around the ring.
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        assert_eq!(
            table_of(&["C"; 6], &bonds),
            vec![1, 2, 1, 3, 1, 4, 2, 5, 3, 6, 4, 5]
        );
    }

    #[test]
    fn ranks_are_a_permutation() {
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        let mut ranks = ranks_of(&["C"; 6], &bonds);
        ranks.sort_unstable();
        assert_eq!(ranks, vec![1, 2, 3, 4, 5, 6]);
    }
}
