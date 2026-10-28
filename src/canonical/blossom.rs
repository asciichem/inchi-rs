//! General graph maximum matching — Edmonds' blossom algorithm.
//! The reference's tautomer decisions are matching-augmentability
//! questions (BalancedNetworkSearch implements Kocay-Stone blossom
//! matching); this is the textbook augmenting-path search,
//! verdict-only — the answer, not the route, is consumed.

#[derive(Debug, Default)]
pub struct Matching {
    adj: Vec<Vec<usize>>,
    matched: Vec<usize>,
}

impl Matching {
    pub fn new(vertices: usize) -> Self {
        Matching {
            adj: vec![Vec::new(); vertices],
            matched: vec![usize::MAX; vertices],
        }
    }

    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.adj[u].push(v);
        self.adj[v].push(u);
    }

    /// Appends a vertex with the given adjacency (the dots).
    pub fn push_vertex(&mut self, neighbors: &[usize]) -> usize {
        let id = self.adj.len();
        self.adj.push(Vec::new());
        self.matched.push(usize::MAX);
        for &n in neighbors {
            self.adj[id].push(n);
            self.adj[n].push(id);
        }
        id
    }

    fn lca(
        matched: &[usize],
        parent: &[usize],
        base: &[usize],
        mut a: usize,
        mut b: usize,
    ) -> usize {
        let mut used = vec![false; matched.len()];
        loop {
            a = base[a];
            used[a] = true;
            if matched[a] == usize::MAX {
                break;
            }
            a = parent[matched[a]];
        }
        loop {
            b = base[b];
            if used[b] {
                return b;
            }
            b = parent[matched[b]];
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn mark_path(
        matched: &[usize],
        parent: &mut [usize],
        base: &[usize],
        blossom: &mut [usize],
        mut v: usize,
        blossom_base: usize,
        mut child: usize,
        marker: usize,
    ) {
        while base[v] != blossom_base {
            blossom[base[v]] = marker;
            blossom[base[matched[v]]] = marker;
            parent[v] = child;
            child = matched[v];
            v = parent[matched[v]];
        }
    }

    /// The augmenting-path search from `root`; applies the flip on
    /// success.
    pub fn augment_from(&mut self, root: usize) -> bool {
        let n = self.adj.len();
        let mut parent = vec![usize::MAX; n];
        let mut base: Vec<usize> = (0..n).collect();
        let mut used = vec![false; n];
        let mut blossom = vec![usize::MAX; n];
        used[root] = true;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(root);
        while let Some(v) = queue.pop_front() {
            for &to in &self.adj[v] {
                if base[v] == base[to] || self.matched[v] == to {
                    continue;
                }
                if to == root || (self.matched[to] != usize::MAX && used[self.matched[to]]) {
                    let curbase = Self::lca(&self.matched, &parent, &base, v, to);
                    for slot in blossom.iter_mut() {
                        *slot = usize::MAX;
                    }
                    Self::mark_path(
                        &self.matched,
                        &mut parent,
                        &base,
                        &mut blossom,
                        v,
                        curbase,
                        to,
                        curbase,
                    );
                    Self::mark_path(
                        &self.matched,
                        &mut parent,
                        &base,
                        &mut blossom,
                        to,
                        curbase,
                        v,
                        curbase,
                    );
                    for vertex in 0..n {
                        if blossom[base[vertex]] == curbase {
                            base[vertex] = curbase;
                            if !used[vertex] {
                                used[vertex] = true;
                                queue.push_back(vertex);
                            }
                        }
                    }
                } else if parent[to] == usize::MAX {
                    parent[to] = v;
                    if self.matched[to] == usize::MAX {
                        let mut w = to;
                        while w != usize::MAX {
                            let pw = parent[w];
                            let nw = self.matched[pw];
                            self.matched[w] = pw;
                            self.matched[pw] = w;
                            w = nw;
                        }
                        return true;
                    } else {
                        used[self.matched[to]] = true;
                        queue.push_back(self.matched[to]);
                    }
                }
            }
        }
        false
    }

    /// A maximum matching; returns the pair count.
    pub fn maximum(&mut self) -> usize {
        for vertex in 0..self.adj.len() {
            if self.matched[vertex] == usize::MAX {
                self.augment_from(vertex);
            }
        }
        self.pairs()
    }

    pub fn pairs(&self) -> usize {
        self.matched.iter().filter(|&&m| m != usize::MAX).count() / 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_triangle_matches_one_pair() {
        let mut m = Matching::new(3);
        m.add_edge(0, 1);
        m.add_edge(1, 2);
        m.add_edge(2, 0);
        assert_eq!(m.maximum(), 1);
    }

    #[test]
    fn a_path_of_four_matches_two_pairs() {
        let mut m = Matching::new(4);
        m.add_edge(0, 1);
        m.add_edge(1, 2);
        m.add_edge(2, 3);
        assert_eq!(m.maximum(), 2);
    }

    #[test]
    fn blossoms_admit_the_odd_cycle_augmentation() {
        // The classic 5-cycle plus pendant: maximum matching 3.
        let mut m = Matching::new(6);
        m.add_edge(0, 1);
        m.add_edge(1, 2);
        m.add_edge(2, 3);
        m.add_edge(3, 4);
        m.add_edge(4, 0);
        m.add_edge(1, 5);
        assert_eq!(m.maximum(), 3);
    }

    #[test]
    fn dot_augmentation_grows_the_matching_once() {
        let mut m = Matching::new(4);
        m.add_edge(0, 1);
        m.add_edge(1, 2);
        m.add_edge(2, 3);
        assert_eq!(m.maximum(), 2);
        let d = m.push_vertex(&[0]);
        let a = m.push_vertex(&[3]);
        assert!(m.augment_from(d));
        assert_eq!(m.pairs(), 3);
        let _ = a;
    }

    #[test]
    fn dots_without_a_route_do_not_augment() {
        let mut m = Matching::new(4);
        m.add_edge(0, 1);
        m.add_edge(2, 3);
        assert_eq!(m.maximum(), 2);
        let d = m.push_vertex(&[0]);
        let _a = m.push_vertex(&[3]);
        assert!(!m.augment_from(d));
    }
}
