//! The linear connection table: for each canonical rank in order, the
//! rank followed by its earlier-ranked neighbors in ascending order.
//! Ordering follows the reference's `CtFullCompare` base layer,
//! including its negative length rule.

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectionTable {
    entries: Vec<u32>,
    atoms: usize,
}

impl ConnectionTable {
    /// Builds the table from final canonical ranks (a permutation of
    /// 1..n).
    pub fn build(ranks: &[u32], adjacency: &[Vec<usize>]) -> ConnectionTable {
        let atoms = adjacency.len();
        let mut entries = Vec::with_capacity(atoms * 3);
        for rank in 1..=atoms as u32 {
            let atom = ranks
                .iter()
                .position(|&r| r == rank)
                .expect("ranks are a permutation of 1..n");
            entries.push(rank);
            let mut earlier: Vec<u32> = adjacency[atom]
                .iter()
                .map(|&b| ranks[b])
                .filter(|&r| r < rank)
                .collect();
            earlier.sort_unstable();
            entries.extend(earlier);
        }
        ConnectionTable { entries, atoms }
    }

    pub fn entries(&self) -> &[u32] {
        &self.entries
    }

    /// The number of atoms the table numbers.
    pub fn atoms(&self) -> usize {
        self.atoms
    }

    /// `CtFullCompare` under the reference's `INCHI_CANON_MIN`
    /// semantics: the search keeps the GREATER result, and greater
    /// means the lexicographically SMALLER connection table. The
    /// reference compares the raw entries, negates the difference
    /// (`INCHI_CANON_MIN` in ichican2.c) and lets the caller replace
    /// its best terminal on a positive result — so a smaller table
    /// wins, and a longer table sorts above a shorter one (the length
    /// difference is negated on entry, then the whole difference is
    /// negated again).
    pub fn compare(&self, other: &ConnectionTable) -> std::cmp::Ordering {
        if self.entries.len() != other.entries.len() {
            return self.entries.len().cmp(&other.entries.len());
        }
        other.entries.cmp(&self.entries)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        adjacency
    }

    #[test]
    fn the_smaller_table_compares_as_greater() {
        // The isobutane skeleton (atom 1 is the hub) under two
        // numberings. Starting at a leaf gives the smaller table.
        let adjacency = graph(4, &[(0, 1), (1, 2), (1, 3)]);
        let leaf_first = ConnectionTable::build(&[2, 1, 3, 4], &adjacency);
        let hub_first = ConnectionTable::build(&[1, 2, 3, 4], &adjacency);
        assert_eq!(leaf_first.entries(), &[1, 2, 1, 3, 1, 4, 1]);
        assert_eq!(hub_first.entries(), &[1, 2, 1, 3, 2, 4, 2]);
        assert_eq!(leaf_first.compare(&hub_first), std::cmp::Ordering::Greater);
        assert_eq!(hub_first.compare(&leaf_first), std::cmp::Ordering::Less);
    }

    #[test]
    fn a_longer_table_compares_as_greater() {
        // The reference negates the length difference on entry and the
        // whole difference again at `done`, so the longer table carries
        // the positive (winning) sign.
        let short = ConnectionTable::build(&[1, 2], &graph(2, &[(0, 1)]));
        let long = ConnectionTable::build(&[1, 2, 3], &graph(3, &[(0, 1), (1, 2)]));
        assert_eq!(long.compare(&short), std::cmp::Ordering::Greater);
        assert_eq!(short.compare(&long), std::cmp::Ordering::Less);
    }
}
