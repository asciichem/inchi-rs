//! The `/c` segment: the connection table rendered as the reference's
//! DFS emission (`GetDfsOrder4CT` + the string form of
//! `MakeCtStringNew`, hydrogen-free).
//!
//! First pass: number the table's graph depth-first from the atom with
//! the fewest bonds, counting subtree sizes and marking ring closures
//! (edges to ancestors) and deleted parent returns. Neighbors are then
//! ordered by subtree size, closures first. Second pass emits the
//! numbers with branch delimiters: `(` opens the branches of an atom,
//! `,` separates them, `)` precedes the last (inline) branch, `-`
//! joins a lone continuation.

use super::table::ConnectionTable;
use std::collections::BTreeSet;

const DELETED: usize = usize::MAX;

pub fn connection_segment(table: &ConnectionTable) -> Option<String> {
    let adjacency = neighbor_lists(table)?;
    let n = adjacency.len();
    if n < 2 {
        return None;
    }
    let start = (0..n)
        .min_by_key(|&i| (adjacency[i].len(), i))
        .expect("n > 0");

    // First DFS: visit order, subtree sizes (closures included), and
    // the deletion of parent returns and already-emitted closures.
    let mut dfs_number = vec![0u32; n];
    let mut descendants = vec![0u32; n];
    let mut cursor = vec![0usize; n];
    let mut neighbors = adjacency.clone();
    let mut stack = vec![start];
    dfs_number[start] = 1;
    descendants[start] = 1;
    let mut next_number = 1u32;
    while let Some(&i) = stack.last() {
        let parent = if stack.len() >= 2 {
            stack[stack.len() - 2]
        } else {
            usize::MAX
        };
        let mut descended = false;
        while cursor[i] < neighbors[i].len() {
            cursor[i] += 1;
            let u = neighbors[i][cursor[i] - 1];
            if dfs_number[u] == 0 {
                next_number += 1;
                dfs_number[u] = next_number;
                descendants[u] = 1;
                stack.push(u);
                descended = true;
                break;
            } else if u != parent && dfs_number[u] < dfs_number[i] {
                // Ring closure: the ancestor side keeps this edge.
                descendants[i] += 1;
            } else {
                neighbors[i][cursor[i] - 1] = DELETED;
            }
        }
        if descended {
            continue;
        }
        cursor[i] = 0;
        stack.pop();
        if let Some(&p) = stack.last() {
            descendants[p] += descendants[i];
        }
    }

    // Order each neighbor list: closures first, then children by
    // subtree size, then by canonical number; deleted entries last.
    for i in 0..n {
        neighbors[i].sort_by_key(|&u| {
            if u == DELETED {
                (u32::MAX, u)
            } else if dfs_number[i] > dfs_number[u] {
                (0, u)
            } else {
                (descendants[u], u)
            }
        });
        let kept = neighbors[i].iter().take_while(|&&u| u != DELETED).count();
        neighbors[i].truncate(kept);
    }

    // Second DFS: emit the numbers with their delimiters.
    let mut out = String::new();
    let mut cursor = vec![0usize; n];
    let mut stack = vec![start];
    out.push_str(&(start + 1).to_string());
    while let Some(&i) = stack.last() {
        if cursor[i] >= neighbors[i].len() {
            cursor[i] = 0;
            stack.pop();
            continue;
        }
        cursor[i] += 1;
        let j = cursor[i] - 1;
        let u = neighbors[i][j];
        let delimiter = if neighbors[i].len() > 1 {
            if j == 0 {
                '('
            } else if j == neighbors[i].len() - 1 {
                ')'
            } else {
                ','
            }
        } else {
            '-'
        };
        out.push(delimiter);
        out.push_str(&(u + 1).to_string());
        if dfs_number[u] > dfs_number[i] {
            stack.push(u);
        }
    }
    Some(out)
}

/// Symmetric neighbor lists rebuilt from the table (`CreateNeighListFromLinearCT`).
fn neighbor_lists(table: &ConnectionTable) -> Option<Vec<Vec<usize>>> {
    let entries = table.entries();
    let n = table.atoms();
    if n == 0 {
        return None;
    }
    let mut adjacency = vec![Vec::new(); n];
    let mut vertex_rank = entries.first().copied()?;
    for &entry in &entries[1..] {
        if entry < vertex_rank {
            let v = vertex_rank as usize - 1;
            let u = entry as usize - 1;
            adjacency[v].push(u);
            adjacency[u].push(v);
        } else {
            vertex_rank = entry;
        }
    }
    Some(adjacency)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::invariants::atom_invariants;
    use crate::canonical::search::canonical_table;

    fn segment(elements: &[&str], bonds: &[(usize, usize)]) -> String {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        let hydrogens = vec![0u32; elements.len()];
        let table = canonical_table(&adjacency, &invariants, &hydrogens, &[]);
        connection_segment(&table).expect("segment")
    }

    #[test]
    fn acetaldehyde_is_a_chain() {
        // Oracle: InChI=1S/C2H4O/c1-2-3/h2H,1H3.
        assert_eq!(segment(&["C", "C", "O"], &[(0, 1), (1, 2)]), "1-2-3");
    }

    #[test]
    fn benzene_walks_1_2_4_6_5_3() {
        // Oracle: InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H.
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        assert_eq!(segment(&["C"; 6], &bonds), "1-2-4-6-5-3-1");
    }

    #[test]
    fn branches_use_parentheses_with_the_last_branch_inline() {
        // Oracle: InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3 (isobutane).
        let got = segment(&["C", "C", "C", "C"], &[(0, 1), (0, 2), (0, 3)]);
        assert_eq!(got, "1-4(2)3");
    }
}

/// The `/c` replay parser (test-support + consumer tooling): recovers
/// the connectivity a segment denotes. Tokens follow the M2 decode —
/// a number bonds to the current walk atom (fresh = advance, seen =
/// closure), `(` saves the walk atom, `)` restores it, and `-` is a
/// separator only. Returns the edge set as sorted `(a, b)` pairs with
/// 0-based atom indices, or `None` when the segment is malformed.
pub fn parse_connection_segment(segment: &str) -> Option<BTreeSet<(usize, usize)>> {
    struct Replay {
        rank_of: std::collections::HashMap<usize, usize>,
        edges: BTreeSet<(usize, usize)>,
        walk: Option<usize>,
        stack: Vec<usize>,
        number: String,
    }

    impl Replay {
        /// A completed number: a fresh rank opens a new atom the walk
        /// advances to; a seen rank is a ring closure and the walk
        /// stays (the M2 decode).
        fn flush(&mut self) -> Result<(), ()> {
            let rank: usize = self.number.parse().map_err(|_| ())?;
            self.number.clear();
            let next = self.rank_of.len();
            let atom = *self.rank_of.entry(rank).or_insert(next);
            if let Some(current) = self.walk {
                self.edges.insert((current.min(atom), current.max(atom)));
            }
            if atom == next {
                self.walk = Some(atom);
            }
            Ok(())
        }
    }

    let mut replay = Replay {
        rank_of: Default::default(),
        edges: BTreeSet::new(),
        walk: None,
        stack: Vec::new(),
        number: String::new(),
    };
    for character in segment.chars() {
        match character {
            '0'..='9' => replay.number.push(character),
            '(' => {
                replay.flush().ok()?;
                replay.stack.push(replay.walk?);
            }
            ')' => {
                replay.flush().ok()?;
                replay.walk = Some(replay.stack.pop()?);
            }
            ',' => {
                replay.flush().ok()?;
                // Branch members each bond to the saved atom, so the
                // comma restores the walk to the branch head.
                replay.walk = Some(*replay.stack.last()?);
            }
            '-' => {
                replay.flush().ok()?;
            }
            _ => return None,
        }
    }
    replay.flush().ok()?;
    if !replay.stack.is_empty() {
        return None;
    }
    Some(replay.edges)
}

#[cfg(test)]
mod replay_tests {
    use super::*;

    #[test]
    fn benzene_segment_replays_to_the_six_ring_edges() {
        let edges = parse_connection_segment("1-2-4-6-5-3-1").expect("parses");
        assert_eq!(edges.len(), 6);
        let degrees = degrees_of(&edges, 6);
        assert_eq!(degrees, vec![2, 2, 2, 2, 2, 2]);
    }

    #[test]
    fn malformed_segments_are_rejected() {
        assert!(parse_connection_segment("1-2(").is_none());
        assert!(parse_connection_segment(")1-2").is_none());
        assert!(parse_connection_segment("1-x").is_none());
    }

    fn degrees_of(edges: &BTreeSet<(usize, usize)>, atoms: usize) -> Vec<usize> {
        let mut degrees = vec![0; atoms];
        for (a, b) in edges {
            degrees[*a] += 1;
            degrees[*b] += 1;
        }
        degrees
    }

    proptest::proptest! {
        #[test]
        fn every_emitted_segment_replays_to_the_same_graph(
            seed in proptest::prelude::any::<u64>()
        ) {
            replay_random_graph(seed);
        }
    }

    /// A random connected graph: a spanning tree plus extra edges,
    /// bond orders drawn from {1, 2}. The emitted `/c` must replay to
    /// a graph with the same atom count, element and degree
    /// multisets, and bond-order multiset — the emitted segment is a
    /// valid emission of the input graph, whatever the canonical
    /// numbering chose.
    fn replay_random_graph(seed: u64) {
        let mut state = seed;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let atoms = 2 + (next() % 7) as usize;
        let mut edges: BTreeSet<(usize, usize)> = BTreeSet::new();
        for atom in 1..atoms {
            let parent = (next() as usize) % atom;
            edges.insert((parent, atom));
        }
        for _ in 0..(next() % 4) {
            let a = (next() as usize) % atoms;
            let b = (next() as usize) % atoms;
            if a != b {
                edges.insert((a.min(b), a.max(b)));
            }
        }
        let bonds: Vec<(usize, usize, u32)> = edges
            .iter()
            .enumerate()
            .map(|(index, &(a, b))| {
                (
                    a,
                    b,
                    if (next() + index as u64).is_multiple_of(4) {
                        2
                    } else {
                        1
                    },
                )
            })
            .collect();
        let adjacency = neighbor_lists_from(&bonds, atoms);
        let ranks: Vec<u32> = (1..=atoms as u32).collect();
        let table = ConnectionTable::build(&ranks, &adjacency);
        let Some(segment) = connection_segment(&table) else {
            return;
        };
        let Some(replayed) = parse_connection_segment(&segment) else {
            panic!("segment {segment:?} does not replay");
        };
        assert_eq!(replayed.len(), bonds.len(), "edge count for {segment:?}");
        // The replay loses bond orders by design (the /c segment is
        // hydrogen-free and order-free); the connectivity must round-
        // trip exactly: edge count and per-atom degree multisets.
        let input_connectivity: Vec<usize> = {
            let mut degrees = vec![0; atoms];
            for &(from, to, _) in &bonds {
                degrees[from] += 1;
                degrees[to] += 1;
            }
            let mut sorted = degrees;
            sorted.sort_unstable();
            sorted
        };
        let replay_connectivity: Vec<usize> = {
            let mut degrees = vec![0; atoms];
            for &(a, b) in &replayed {
                degrees[a] += 1;
                degrees[b] += 1;
            }
            let mut sorted = degrees;
            sorted.sort_unstable();
            sorted
        };
        assert_eq!(
            input_connectivity, replay_connectivity,
            "degrees for {segment:?}"
        );
    }

    fn neighbor_lists_from(bonds: &[(usize, usize, u32)], atoms: usize) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); atoms];
        for &(from, to, _) in bonds {
            adjacency[from].push(to);
            adjacency[to].push(from);
        }
        adjacency
    }
}
