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

    /// `CtFullCompare` layer 0: the first difference decides, and a
    /// longer table compares as the smaller one.
    pub fn compare(&self, other: &ConnectionTable) -> std::cmp::Ordering {
        if self.entries.len() != other.entries.len() {
            return other.entries.len().cmp(&self.entries.len());
        }
        self.entries.cmp(&other.entries)
    }
}
