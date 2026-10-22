//! Initial atom invariants: the reference's element ordering (carbon
//! first, then the other elements by their two-byte symbol key, H
//! last) and the heavy-atom degree. Deliberately hydrogen-free: the
//! oracle's `/c` numbering is decided before hydrogen counts enter.

/// Sort key standing in for the reference's two-byte element field:
/// one-byte symbols pad with a space, matching its byte order.
fn element_key(element: &str) -> [u8; 2] {
    let bytes = element.as_bytes();
    [
        bytes.first().copied().unwrap_or(b' '),
        bytes.get(1).copied().unwrap_or(b' '),
    ]
}

pub fn atom_invariants(elements: &[String], adjacency: &[Vec<usize>]) -> Vec<(u16, u16)> {
    let mut others: Vec<[u8; 2]> = elements
        .iter()
        .filter(|e| e.as_str() != "C" && e.as_str() != "H")
        .map(|e| element_key(e))
        .collect();
    others.sort_unstable();
    others.dedup();
    let rank = |element: &str| -> u16 {
        if element == "C" {
            1
        } else {
            2 + others
                .iter()
                .position(|key| *key == element_key(element))
                .unwrap_or_else(|| {
                    others
                        .iter()
                        .position(|key| *key > element_key(element))
                        .unwrap_or(others.len())
                }) as u16
        }
    };
    elements
        .iter()
        .enumerate()
        .map(|(index, element)| (rank(element), adjacency[index].len() as u16))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn carbon_sorts_before_the_other_elements() {
        let elements: Vec<String> = ["N", "C", "Cl", "Br", "B", "Si"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let adjacency: Vec<Vec<usize>> = elements.iter().map(|_| Vec::new()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        // C first; then B, Br, Cl, N, Si by byte order.
        let ranks: Vec<u16> = invariants.iter().map(|(r, _)| *r).collect();
        assert_eq!(ranks, vec![5, 1, 4, 3, 2, 6]);
    }
}
