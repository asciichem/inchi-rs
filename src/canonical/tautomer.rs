//! Mobile-H (tautomer) detection — M3, ported from the reference's
//! ichitaut semantics. An endpoint is O/S/Se/Te (endpoint valence 2)
//! or N (3) whose attached hydrogens plus bond-order sum exactly fill
//! the endpoint valence; the fill decides the role — no multiple bond
//! = donor (the H's owner), exactly one double-bond worth = acceptor.
//! Overflowing the endpoint valence disqualifies an atom outright, so
//! a nitro nitrogen (bond sum 4 over N's 3) keeps its aci hydrogens
//! fixed and a sulfone/sulfate sulfur (bond sum ≥ 3 over S's 2)
//! never participates — only thiocarbonyl =S and S-H do.
//! Endpoints link through a shared anchor that holds both a donor and
//! an acceptor (the 1,3 template, the donor's bond to the anchor a
//! real single bond carrying its H), and along alternating
//! single/aromatic/double shift walks whose far terminal is the
//! COMPLEMENTARY role — a phenolic O can never reach another phenolic
//! O across an all-carbon aromatic arc because both terminals are
//! donors. An aromatic run inside a walk is entered from the donor's
//! own bond or a double bond and left into a double bond or the donor
//! terminal: shift paths never jump π systems through an interior
//! chain single. The canonicalization ranks atoms by hydrogens with mobile
//! H removed (the reference's numHNoTautH — see
//! [`without_mobile_h`]), which is what settles terminal-O ordering;
//! a group's hydrogens are the ±1-standardized counts of its members.

use crate::graph::Structure;
use std::collections::{BTreeSet, HashMap};

/// One mobile-H group: the hydrogen count and the canonical ranks of
/// its endpoints (ascending).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MobileGroup {
    pub count: u32,
    pub ranks: Vec<u32>,
}

/// Endpoint valence per the reference's `get_endpoint_valence`.
fn endpoint_valence(element: &str) -> Option<u32> {
    match element {
        "O" | "S" | "Se" | "Te" => Some(2),
        "N" => Some(3),
        _ => None,
    }
}

/// How an endpoint's bond-order sum relates to its heavy-atom valence:
/// all single bonds = [`Role::Donor`], one double-bond worth =
/// [`Role::Acceptor`].
#[derive(Clone, Copy, PartialEq, Eq)]
enum Role {
    Donor,
    Acceptor,
}

/// Cached per-atom structure facts the endpoint arithmetic needs.
struct Context<'a> {
    structure: &'a Structure,
    hydrogens: &'a [u32],
    neighbors: Vec<Vec<(usize, u32)>>,
    bond_index: HashMap<(usize, usize), usize>,
    aromatic: BTreeSet<usize>,
    /// Bond-order sums, doubled so aromatic 1.5 stays integral.
    bond_sum2: Vec<i32>,
}

impl<'a> Context<'a> {
    fn new(structure: &'a Structure, hydrogens: &'a [u32]) -> Self {
        let n = structure.atoms.len();
        let aromatic = crate::canonical::aromatic::aromatic_bonds(structure);
        let mut neighbors = vec![Vec::new(); n];
        let mut bond_index = HashMap::new();
        let mut bond_sum2 = vec![0i32; n];
        for (index, &(from, to, kind)) in structure.bonds.iter().enumerate() {
            let order2 = if aromatic.contains(&index) {
                3
            } else {
                match kind {
                    2 => 4,
                    3 => 6,
                    _ => 2,
                }
            };
            neighbors[from].push((to, kind));
            neighbors[to].push((from, kind));
            bond_index.insert((from.min(to), from.max(to)), index);
            bond_sum2[from] += order2;
            bond_sum2[to] += order2;
        }
        Context {
            structure,
            hydrogens,
            neighbors,
            bond_index,
            aromatic,
            bond_sum2,
        }
    }

    fn is_aromatic(&self, a: usize, b: usize) -> bool {
        self.bond_index
            .get(&(a.min(b), a.max(b)))
            .map(|&index| self.aromatic.contains(&index))
            .unwrap_or(false)
    }

    /// The reference's `nGetEndpointInfo`.
    fn endpoint(&self, atom: usize) -> Option<Role> {
        let record = &self.structure.atoms[atom];
        let valence = endpoint_valence(&record.element)? as i32;
        let heavy = self.neighbors[atom].len() as i32;
        if valence <= heavy {
            return None;
        }
        if record.charge > 0 {
            return None;
        }
        let mobile = self.mobile(atom) as i32;
        if 2 * mobile + self.bond_sum2[atom] != 2 * valence {
            return None;
        }
        match self.bond_sum2[atom] - 2 * heavy {
            0 => Some(Role::Donor),
            2 => Some(Role::Acceptor),
            _ => None,
        }
    }

    /// Mobile H the reference's `nMobile`. The ±1-standardized counts
    /// already carry the anion's riding proton (an anionic O shows H 1
    /// exactly where the reference shows num_H 0 plus one negative
    /// charge), so the stored count is the whole story.
    fn mobile(&self, atom: usize) -> u32 {
        self.hydrogens[atom]
    }
}

/// The atoms whose hydrogens are mobile: members of a mobile group.
/// Endpoint arithmetic alone is not enough — an alcohol oxygen is a
/// donor candidate that never joins a group, and its hydrogen must
/// keep counting in the canonicalization.
pub fn endpoint_mask(structure: &Structure, hydrogens: &[u32]) -> Vec<bool> {
    let context = Context::new(structure, hydrogens);
    let mut mask = vec![false; structure.atoms.len()];
    for members in components(&context, hydrogens) {
        for atom in members {
            mask[atom] = true;
        }
    }
    mask
}

/// Stage-2 hydrogen counts: the given counts with mobile H removed —
/// the reference's `numHNoTautH`, which both the `/c` and `/h`
/// canonicalizations rank by.
pub fn without_mobile_h(hydrogens: &[u32], mask: &[bool]) -> Vec<u32> {
    hydrogens
        .iter()
        .zip(mask)
        .map(|(&count, &endpoint)| if endpoint { 0 } else { count })
        .collect()
}

fn build_links(context: &Context) -> Vec<BTreeSet<usize>> {
    let n = context.structure.atoms.len();
    let mut linked: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];

    // The 1,3 template: an anchor holding both a donor and an
    // acceptor joins all of its endpoint neighbors into one group.
    // A real single bond counts toward the donor role only when the
    // endpoint actually carries mobile H; an aromatic bond takes the
    // endpoint's arithmetic role as either.
    for anchor in 0..n {
        let mut ends: Vec<usize> = Vec::new();
        let mut donors = 0;
        let mut acceptors = 0;
        let mut any_mobile = false;
        for &(neighbor, kind) in &context.neighbors[anchor] {
            let Some(role) = context.endpoint(neighbor) else {
                continue;
            };
            let aromatic = context.is_aromatic(anchor, neighbor);
            let classified = if aromatic {
                role
            } else {
                match kind {
                    1 if context.mobile(neighbor) > 0 => Role::Donor,
                    2 => Role::Acceptor,
                    _ => continue,
                }
            };
            match classified {
                Role::Donor => donors += 1,
                Role::Acceptor => acceptors += 1,
            }
            any_mobile |= context.mobile(neighbor) > 0;
            ends.push(neighbor);
        }
        if ends.len() > 1 && donors > 0 && acceptors > 0 && any_mobile {
            for &x in &ends {
                for &y in &ends {
                    if x != y {
                        linked[x].insert(y);
                        linked[y].insert(x);
                    }
                }
            }
        }
    }

    // Shift walks from each endpoint toward its complementary role:
    // a donor walks single-bond-first (the H leaves along its own
    // single bond), an acceptor walks double/aromatic-bond-first.
    // The 1,5 template terminates on a double bond — a bare
    // single-terminated acceptor is not one (phenol OH's stay fixed).
    for start in 0..n {
        let Some(role) = context.endpoint(start) else {
            continue;
        };
        match role {
            Role::Donor => {
                walk_path(
                    context,
                    start,
                    start,
                    1,
                    Far::Acceptor,
                    2,
                    None,
                    &mut linked,
                );
                walk_path(
                    context,
                    start,
                    start,
                    1,
                    Far::Acceptor,
                    4,
                    None,
                    &mut linked,
                );
            }
            Role::Acceptor => {
                walk_path(context, start, start, 2, Far::Donor, 2, None, &mut linked);
                walk_path(context, start, start, 2, Far::Donor, 4, None, &mut linked);
            }
        }
    }
    linked
}

/// What the far terminal of a shift walk must be: the complementary
/// endpoint role, reached along a bond of its own side (double-like
/// into the acceptor, a real single into the donor).
#[derive(Clone, Copy, PartialEq)]
enum Far {
    Acceptor,
    Donor,
}

/// The π class of a walked bond: an aromatic bond delocalizes the
/// shift, a real double carries it, a real single is the donor's own
/// terminal or the chain link between π centers.
#[derive(Clone, Copy, PartialEq)]
enum BondClass {
    Aromatic,
    Single,
    Double,
}

/// A real single bonds to an aromatic run only as the donor's own
/// first bond or the far donor's terminal — the reference's shift
/// paths never jump off a π system through an interior chain single
/// (the aryl-ester carbonyl that reaches across its ring is the
/// canonical reject).
fn transition_leaves_pi_system(previous: &Step, current: BondClass, last_bond: bool) -> bool {
    match previous {
        // The walk origin's own terminal bond is exempt — the donor
        // leaves along its single into the π system.
        Step::OwnTerminal => false,
        Step::Over(BondClass::Single) => current == BondClass::Aromatic,
        Step::Over(BondClass::Aromatic) => current == BondClass::Single && !last_bond,
        Step::Over(_) => false,
    }
}

/// The bond just walked: the origin's own terminal, or an interior
/// bond of the chain.
#[derive(Clone, Copy)]
enum Step {
    OwnTerminal,
    Over(BondClass),
}

#[allow(clippy::too_many_arguments)]
fn walk_path(
    context: &Context,
    atom: usize,
    start: usize,
    parity: u32,
    far: Far,
    remaining: usize,
    previous: Option<Step>,
    linked: &mut Vec<BTreeSet<usize>>,
) {
    if remaining == 0 {
        return;
    }
    for (next, kind) in &context.neighbors[atom] {
        let next = *next;
        let aromatic_bond = context.is_aromatic(atom, next);
        if *kind != parity && !aromatic_bond {
            continue;
        }
        let class = if aromatic_bond {
            BondClass::Aromatic
        } else if *kind == 1 {
            BondClass::Single
        } else {
            BondClass::Double
        };
        if let Some(previous) = &previous {
            if transition_leaves_pi_system(previous, class, remaining == 1) {
                continue;
            }
        }
        if remaining == 1 {
            if next == start || context.endpoint(next) != required_role(far) {
                continue;
            }
            let reached = match far {
                // The walk's far acceptor is met by a real double — a
                // carbonyl met across an aromatic run is not a shift
                // target (the exocyclic quinone C=O is; a ring N
                // reached through the ring pairs only via a shared
                // anchor).
                Far::Acceptor => class == BondClass::Double,
                Far::Donor => class == BondClass::Single,
            };
            if reached {
                linked[start].insert(next);
                linked[next].insert(start);
            }
            continue;
        }
        let step = match previous {
            None => Step::OwnTerminal,
            Some(_) => Step::Over(class),
        };
        walk_path(
            context,
            next,
            start,
            if parity == 1 { 2 } else { 1 },
            far,
            remaining - 1,
            Some(step),
            linked,
        );
    }
}

fn required_role(far: Far) -> Option<Role> {
    match far {
        Far::Acceptor => Some(Role::Acceptor),
        Far::Donor => Some(Role::Donor),
    }
}

/// Groups of endpoints connected by 1,3 anchors or shift walks,
/// filtered to those with at least two members and at least one
/// mobile H. Rank-free — membership depends only on the structure.
fn components(context: &Context, hydrogens: &[u32]) -> Vec<Vec<usize>> {
    let n = context.structure.atoms.len();
    let linked = build_links(context);
    let mut component = vec![usize::MAX; n];
    let mut groups = Vec::new();
    for start in 0..n {
        if context.endpoint(start).is_none() || component[start] != usize::MAX {
            continue;
        }
        let mut members = Vec::new();
        let mut queue = vec![start];
        component[start] = start;
        while let Some(atom) = queue.pop() {
            members.push(atom);
            for &next in &linked[atom] {
                if component[next] == usize::MAX {
                    component[next] = start;
                    queue.push(next);
                }
            }
        }
        let count: u32 = members.iter().map(|&atom| hydrogens[atom]).sum();
        if members.len() >= 2 && count > 0 {
            groups.push(members);
        }
    }
    groups
}

pub fn mobile_groups(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<MobileGroup> {
    let context = Context::new(structure, hydrogens);
    let mut groups = Vec::new();
    for members in components(&context, hydrogens) {
        let count: u32 = members.iter().map(|&atom| hydrogens[atom]).sum();
        let mut group_ranks: Vec<u32> = members.iter().map(|&atom| ranks[atom]).collect();
        group_ranks.sort_unstable();
        groups.push(MobileGroup {
            count,
            ranks: group_ranks,
        });
    }
    groups.sort_by_key(|group| (group.ranks.len(), group.ranks[0]));
    groups
}

/// The atoms whose hydrogens live in mobile groups (excluded from the
/// fixed-H entries).
pub fn endpoint_atoms(ranks: &[u32], groups: &[MobileGroup]) -> BTreeSet<usize> {
    let mut by_rank = std::collections::BTreeMap::new();
    for (atom, &rank) in ranks.iter().enumerate() {
        by_rank.insert(rank, atom);
    }
    let mut endpoints = BTreeSet::new();
    for group in groups {
        for &rank in &group.ranks {
            if let Some(&atom) = by_rank.get(&rank) {
                endpoints.insert(atom);
            }
        }
    }
    endpoints
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
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure { atoms, bonds }
    }

    #[test]
    fn carboxyl_pairs_the_two_oxygens() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
    }

    #[test]
    fn alcohol_oxygen_is_not_an_endpoint() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O")],
            vec![(0, 1, 1), (1, 2, 1)],
        );
        assert!(endpoint_mask(&s, &[3, 2, 1]).iter().all(|&e| !e));
        assert!(mobile_groups(&s, &[3, 2, 1], &[1, 2, 3]).is_empty());
    }

    #[test]
    fn amide_nitrogen_joins_the_carbonyl_oxygen() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("N")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 2], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 2);
        assert_eq!(groups[0].ranks, vec![3, 4]);
        let endpoints = endpoint_atoms(&[1, 2, 3, 4], &groups);
        assert_eq!(endpoints, BTreeSet::from([2, 3]));
    }

    #[test]
    fn the_nitro_anion_groups_its_oxygens_but_not_its_nitrogen() {
        // CH3-N(=O)-O(-): the oracle records (H,3,4) over the two
        // oxygens; the nitrogen's bond sum (4) overflows the endpoint
        // valence of N (3), so its hydrogens stay fixed.
        let s = structure(
            vec![atom("C"), atom("N"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[3, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
        assert!(!endpoint_mask(&s, &[3, 0, 0, 1])[1]);
    }

    #[test]
    fn sulfone_sulfur_never_participates_but_thiocarbonyl_does() {
        // R-S(=O)2-R': bond sum 6 overflows S's endpoint valence 2.
        let s = structure(
            vec![atom("C"), atom("S"), atom("O"), atom("O"), atom("C")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 2), (1, 4, 1)],
        );
        assert!(!endpoint_mask(&s, &[3, 0, 0, 0, 3])[1]);
        // R-C(=S)-NH2: the thiocarbonyl S is an acceptor endpoint.
        let s = structure(
            vec![atom("C"), atom("S"), atom("N")],
            vec![(0, 1, 2), (0, 2, 1)],
        );
        assert!(endpoint_mask(&s, &[0, 0, 2])[1]);
    }

    #[test]
    fn aromatic_bonds_carry_the_shift_across_a_ring_to_a_ring_acceptor() {
        // 2-hydroxypyridine skeleton: the OH oxygen (donor) and the
        // ring N (acceptor) are tautomeric through the aromatic ring.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms[0] = atom("N");
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 1), // ring N - C1
                (1, 2, 2),
                (2, 3, 1),
                (3, 4, 2),
                (4, 5, 1),
                (5, 0, 2), // ring closes with a double to N
                (1, 6, 1), // exocyclic OH on C1
            ],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 0, 0, 0, 1], &[1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![1, 7]);
    }

    #[test]
    fn two_phenolic_oxygens_never_merge_across_an_all_carbon_ring() {
        // Catechol skeleton: both O are donors, and a shift walk
        // needs complementary terminals — no group.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms.push(atom("O"));
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 1),
                (1, 2, 2),
                (2, 3, 1),
                (3, 4, 2),
                (4, 5, 1),
                (5, 0, 2),
                (1, 6, 1),
                (2, 7, 1),
            ],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 0, 0, 0, 1, 1], &[1, 2, 3, 4, 5, 6, 7, 8]);
        assert!(groups.is_empty());
    }
}
