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

/// The ring perception of the reference's MarkRingSystemsInp: a DFS
/// segmentation at `dfs == low` gives every atom a ring-system id and
/// size — a chain atom is a size-1 system of its own, rings joined by
/// a bridge are separate systems, spiro fusions are one — plus the
/// articulation points the 1,4 template requires.
struct RingSystems {
    system: Vec<usize>,
    size: Vec<usize>,
    cut_vertex: Vec<bool>,
}

fn ring_systems(structure: &Structure) -> RingSystems {
    let n = structure.atoms.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    let mut system = vec![0usize; n];
    let mut size = vec![0usize; n];
    let mut cut_vertex = vec![false; n];

    // Articulation points: a non-root atom whose subtree cannot reach
    // above it; the root only when it has two independent children.
    {
        let mut disc = vec![0u32; n];
        #[allow(clippy::too_many_arguments)]
        fn visit(
            atom: usize,
            parent: Option<usize>,
            parent_edge_used: &mut bool,
            neighbors: &[Vec<usize>],
            disc: &mut Vec<u32>,
            low: &mut Vec<u32>,
            cut_vertex: &mut Vec<bool>,
            timer: &mut u32,
            root: usize,
            root_children: &mut u32,
        ) -> u32 {
            *timer += 1;
            disc[atom] = *timer;
            let mut low_here = *timer;
            for &next in &neighbors[atom] {
                if Some(next) == parent && !*parent_edge_used {
                    *parent_edge_used = true;
                    continue;
                }
                if disc[next] == 0 {
                    *root_children += (atom == root) as u32;
                    let child_low = visit(
                        next,
                        Some(atom),
                        &mut false,
                        neighbors,
                        disc,
                        low,
                        cut_vertex,
                        timer,
                        root,
                        root_children,
                    );
                    low_here = low_here.min(child_low);
                    if child_low >= disc[atom] && (atom != root || *root_children > 1) {
                        cut_vertex[atom] = true;
                    }
                } else {
                    low_here = low_here.min(disc[next]);
                }
            }
            low[atom] = low_here;
            low_here
        }
        for root in 0..n {
            if disc[root] == 0 {
                let mut root_children = 0u32;
                visit(
                    root,
                    None,
                    &mut false,
                    &neighbors,
                    &mut disc,
                    &mut vec![0; n],
                    &mut cut_vertex,
                    &mut 0,
                    root,
                    &mut root_children,
                );
            }
        }
    }

    // Ring-system segmentation: on retracting from an atom whose low
    // equals its own dfs, the DFS stack suffix from it unwinds into a
    // numbered system of that size.
    {
        let mut disc = vec![0u32; n];
        let mut timer = 0u32;
        let mut stack: Vec<usize> = Vec::new();
        let mut next_system = 0usize;
        #[allow(clippy::too_many_arguments)]
        fn visit(
            atom: usize,
            parent: Option<usize>,
            parent_edge_used: &mut bool,
            neighbors: &[Vec<usize>],
            disc: &mut Vec<u32>,
            timer: &mut u32,
            stack: &mut Vec<usize>,
            system: &mut Vec<usize>,
            size: &mut Vec<usize>,
            next_system: &mut usize,
        ) -> u32 {
            *timer += 1;
            disc[atom] = *timer;
            let mut low_here = *timer;
            stack.push(atom);
            for &next in &neighbors[atom] {
                if Some(next) == parent && !*parent_edge_used {
                    *parent_edge_used = true;
                    continue;
                }
                let child_low = if disc[next] == 0 {
                    visit(
                        next,
                        Some(atom),
                        &mut false,
                        neighbors,
                        disc,
                        timer,
                        stack,
                        system,
                        size,
                        next_system,
                    )
                } else {
                    disc[next]
                };
                low_here = low_here.min(child_low);
            }
            if disc[atom] == low_here {
                *next_system += 1;
                let mut count = 0usize;
                while let Some(top) = stack.pop() {
                    system[top] = *next_system;
                    count += 1;
                    if top == atom {
                        break;
                    }
                }
                size[atom] = count;
                let top = system.len();
                for index in 0..top {
                    if system[index] == *next_system {
                        size[index] = count;
                    }
                }
            }
            low_here
        }
        for root in 0..n {
            if disc[root] == 0 {
                visit(
                    root,
                    None,
                    &mut false,
                    &neighbors,
                    &mut disc,
                    &mut timer,
                    &mut stack,
                    &mut system,
                    &mut size,
                    &mut next_system,
                );
            }
        }
    }

    RingSystems {
        system,
        size,
        cut_vertex,
    }
}

/// Cached per-atom structure facts the endpoint arithmetic needs.
struct Context<'a> {
    structure: &'a Structure,
    hydrogens: &'a [u32],
    neighbors: Vec<Vec<(usize, u32)>>,
    bond_index: HashMap<(usize, usize), usize>,
    /// Perceived-aromatic bonds: the reference's normalization marks
    /// them BOND_ALTERN while the valence arithmetic keeps the input
    /// kekule's integer orders.
    aromatic: BTreeSet<usize>,
    /// Bond-order sums of the INPUT kekule orders, doubled.
    bond_sum2: Vec<i32>,
    /// The MarkRingSystemsInp perception the templates scan.
    rings: RingSystems,
}

impl<'a> Context<'a> {
    fn new(structure: &'a Structure, hydrogens: &'a [u32]) -> Self {
        let n = structure.atoms.len();
        let aromatic = crate::canonical::aromatic::aromatic_bonds(structure);
        let mut neighbors = vec![Vec::new(); n];
        let mut bond_index = HashMap::new();
        let mut bond_sum2 = vec![0i32; n];
        for (index, &(from, to, kind)) in structure.bonds.iter().enumerate() {
            // The INPUT kekule orders carry the valence arithmetic.
            let order2 = match kind {
                2 => 4,
                3 => 6,
                4 => 3,
                _ => 2,
            };
            neighbors[from].push((to, kind));
            neighbors[to].push((from, kind));
            bond_index.insert((from.min(to), from.max(to)), index);
            bond_sum2[from] += order2;
            bond_sum2[to] += order2;
        }
        let rings = ring_systems(structure);
        Context {
            structure,
            hydrogens,
            neighbors,
            bond_index,
            aromatic,
            bond_sum2,
            rings,
        }
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

    /// The reference reads `bond_type & ~BOND_MARK_ALL`, and
    /// SetTautomericBonds REPLACES the marked bond's type with
    /// BOND_TAUTOM — a registered bond reads TAUTOM, its original
    /// order gone. A perceived-aromatic bond reads BOND_ALTERN.
    fn bond_type(&self, state: &Registration, a: usize, b: usize) -> BondType {
        let &index = self
            .bond_index
            .get(&(a.min(b), a.max(b)))
            .expect("template walk on a bonded pair");
        if state.marks.contains(&index) {
            return BondType::Tautom;
        }
        if self.aromatic.contains(&index) {
            return BondType::Altern;
        }
        match self.structure.bonds[index].2 {
            2 => BondType::Double,
            3 => BondType::Triple,
            4 => BondType::Altern,
            _ => BondType::Single,
        }
    }

    /// The BNS traversal rule the alt path's intermediates answer to:
    /// a relay atom needs valence slack — its bond orders plus H must
    /// sit below its standard valence — or the augmentation cannot
    /// pass (a saturated sp3 CH2 blocks a 1,5 shift; an aromatic or
    /// sp2 atom relays it).
    fn has_valence_slack(&self, atom: usize) -> bool {
        let max = match self.structure.atoms[atom].element.as_str() {
            "C" => 4,
            "N" => 3,
            "O" => 2,
            "P" => 5,
            "S" => 6,
            _ => return false,
        };
        2 * max > self.bond_sum2[atom] + 2 * self.hydrogens[atom] as i32
    }

    /// The reference's `bIsCenterPointStrict`: all-single with an
    /// endpoint valence above it and a mobile H (or an anionic site),
    /// or exactly one double-bond worth on a strict center element.
    /// The reference's `c_point` disjunct is a charged-BNS context
    /// flag our scope never sets.
    fn is_centerpoint_strict(&self, atom: usize) -> bool {
        let record = &self.structure.atoms[atom];
        let heavy = self.neighbors[atom].len() as i32;
        if self.bond_sum2[atom] == 2 * heavy {
            let ev = endpoint_valence(&record.element).unwrap_or(0) as i32;
            return ev != 0 && ev > heavy && (record.h_count > 0 || record.charge == -1);
        }
        if self.bond_sum2[atom] == 2 * (heavy + 1) {
            return matches!(record.element.as_str(), "C" | "N" | "P" | "As" | "Sb");
        }
        false
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

/// The mobile groups rank-free, as (hydrogen count, member atoms) —
/// the terminal-compare rows inside the canonicalization search (the
/// reference's t-group CT section).
pub fn group_membership(structure: &Structure, hydrogens: &[u32]) -> Vec<(u32, Vec<usize>)> {
    let context = Context::new(structure, hydrogens);
    components(&context, hydrogens)
        .into_iter()
        .map(|members| {
            let count = members.iter().map(|&m| hydrogens[m]).sum();
            (count, members)
        })
        .collect()
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

/// The registration state the templates mutate: group membership per
/// atom (the reference's `at[].endpoint`) and the bonds marked
/// tautomeric (`BOND_TAUTOM`). The marks are the transitive device —
/// registering a group marks all its anchors' bonds, and a
/// tautomeric bond lets its endpoint qualify as either role in
/// later scans.
#[derive(Default)]
struct Registration {
    group_of: Vec<Option<usize>>,
    marks: BTreeSet<usize>,
}

impl Registration {
    fn registered(&self, atom: usize) -> bool {
        self.group_of[atom].is_some()
    }

    /// Merges the members into one group; `true` when the state
    /// changed (a member was unregistered or two groups fused).
    fn register(&mut self, members: &[usize]) -> bool {
        let existing: Vec<usize> = members.iter().filter_map(|&m| self.group_of[m]).collect();
        match existing.len() {
            0 => {
                let id = self.group_of.iter().flatten().max().map_or(0, |m| m + 1);
                for &m in members {
                    self.group_of[m] = Some(id);
                }
                true
            }
            1 => {
                let id = existing[0];
                let mut changed = false;
                for &m in members {
                    if self.group_of[m] != Some(id) {
                        self.group_of[m] = Some(id);
                        changed = true;
                    }
                }
                changed
            }
            _ => {
                let keep = *existing.iter().min().unwrap();
                let absorbed: BTreeSet<usize> = existing.into_iter().collect();
                let mut changed = false;
                // Absorb every member of the fused groups, then add
                // the new members themselves.
                for group in self.group_of.iter_mut() {
                    if let Some(id) = *group {
                        if absorbed.contains(&id) && id != keep {
                            *group = Some(keep);
                            changed = true;
                        }
                    }
                }
                for &m in members {
                    if self.group_of[m] != Some(keep) {
                        self.group_of[m] = Some(keep);
                        changed = true;
                    }
                }
                changed
            }
        }
    }
}

/// Builds the endpoint linkage: the templates run as registration
/// rounds — anchors, then the shift walks — each registering groups
/// and re-marking bonds, to a fixpoint. The returned adjacency joins
/// every pair of co-members.
fn build_links(context: &Context) -> Vec<BTreeSet<usize>> {
    let n = context.structure.atoms.len();
    let mut state = Registration {
        group_of: vec![None; n],
        marks: BTreeSet::new(),
    };
    for _ in 0..8 {
        let mut changed = registration_round(context, &mut state);
        if !changed {
            changed = walk_round(context, &mut state);
        }
        if !changed {
            break;
        }
    }
    // Adjacency over co-membership (components() then yields the
    // groups; empty ones drop).
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (atom, group) in state.group_of.iter().enumerate() {
        if let Some(group) = group {
            members[*group].push(atom);
        }
    }
    let mut linked: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for group in members {
        for &x in &group {
            for &y in &group {
                if x != y {
                    linked[x].insert(y);
                    linked[y].insert(x);
                }
            }
        }
    }
    linked
}

/// One registration round: the 1,3 anchors (with the tautomeric
/// marks feeding the classification) and the shift walks.
fn registration_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut changed = false;

    // The 1,3 template: an endpoint reaches its centerpoint through
    // a π bond (double, aromatic, or already-marked tautomeric); the
    // centerpoint's other endpoint neighbors join. A single bond is
    // the donor side only with mobile H or membership; a double bond
    // is the acceptor side; a π/tautomeric bond reads either role
    // from the arithmetic or from membership. Registering marks ALL
    // the centerpoint's bonds — the transitive device that lets the
    // next round admit the mark's other end.
    for anchor in 0..n {
        if context.structure.atoms[anchor].charge != 0 {
            continue;
        }
        // The anchor scans only when some endpoint reaches it
        // through a π-class bond.
        let mut reachable = false;
        for &(neighbor, kind) in &context.neighbors[anchor] {
            if context.endpoint(neighbor).is_some()
                && (kind == 2
                    || context.is_aromatic(anchor, neighbor)
                    || context
                        .bond_index
                        .get(&(anchor.min(neighbor), anchor.max(neighbor)))
                        .is_some_and(|bond| state.marks.contains(bond)))
            {
                reachable = true;
                break;
            }
        }
        if !reachable {
            continue;
        }
        let mut ends: Vec<usize> = Vec::new();
        let mut donors = 0;
        let mut acceptors = 0;
        let mut any_mobile = false;
        for &(neighbor, kind) in &context.neighbors[anchor] {
            let Some(arithmetic) = context.endpoint(neighbor) else {
                continue;
            };
            let registered = state.registered(neighbor);
            let marked = state
                .marks
                .contains(&context.bond_index[&(anchor.min(neighbor), anchor.max(neighbor))]);
            let aromatic = context.is_aromatic(anchor, neighbor);
            let (classified, possible) = if aromatic || marked {
                // An ALTERN bond (perceived aromatic, or a bond the
                // previous round marked) carries either role through
                // the arithmetic itself: donor iff cDonor, acceptor
                // iff cNeutralBondsValence > valence.
                (arithmetic, true)
            } else {
                match kind {
                    1 => {
                        let donor = context.mobile(neighbor) > 0 || registered;
                        (Role::Donor, donor)
                    }
                    2 => (Role::Acceptor, true),
                    _ => continue,
                }
            };
            if !possible {
                continue;
            }
            match classified {
                Role::Donor => donors += 1,
                Role::Acceptor => acceptors += 1,
            }
            any_mobile |= context.mobile(neighbor) > 0 || registered;
            ends.push(neighbor);
        }
        if ends.len() > 1 && donors > 0 && acceptors > 0 && any_mobile {
            changed |= state.register(&ends);
            // SetTautomericBonds marks only the bonds to admitted
            // endpoints, never the centerpoint's other bonds.
            for &end in &ends {
                state
                    .marks
                    .insert(context.bond_index[&(anchor.min(end), anchor.max(end))]);
            }
        }
    }

    changed
}

/// A template check: receives the live registration state, the path
/// atoms, and the path bonds (bond i joins atom i to atom i + 1; the
/// last entry of a ring path is the closure bond).
type DfsCheck<'ck> = &'ck mut dyn FnMut(&mut Registration, &[usize], &[BondType]) -> bool;

/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the walks and ring templates). From each
/// endpoint toward its complementary role: a donor walks
/// single-bond-first (the H leaves along its own single bond), an
/// acceptor walks double/aromatic-bond-first. The 1,5 template
/// terminates on a double bond — a bare single-terminated acceptor
/// is not one (phenol OH's stay fixed) — and needs one endpoint
/// OFF-ring (TAUT_15_NON_RING): two ring-bound endpoints pair only
/// through their own ring templates.
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the TAUT_OTHER block). Four template
/// families, in the reference's order: the 1,5 alt path
/// (TAUT_15_NON_RING), 6-member alt rings (TAUT_4PYRIDINOL_RINGS),
/// pyrazole 5-rings (TAUT_PYRAZOLE_RINGS), and the 1,4 tropolone
/// pattern (TAUT_TROPOLONE_5/7).
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the walks and ring templates). From each
/// endpoint toward its complementary role: a donor walks
/// single-bond-first (the H leaves along its own single bond), an
/// acceptor walks double/aromatic-bond-first. The 1,5 template
/// terminates on a double bond — a bare single-terminated acceptor
/// is not one (phenol OH's stay fixed) — and needs one endpoint
/// OFF-ring (TAUT_15_NON_RING): two ring-bound endpoints pair only
/// through their own ring templates.
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the TAUT_OTHER block). Four template
/// families, in the reference's order: the 1,5 alt path
/// (TAUT_15_NON_RING), 6-member alt rings (TAUT_4PYRIDINOL_RINGS),
/// pyrazole 5-rings (TAUT_PYRAZOLE_RINGS), and the 1,4 tropolone
/// pattern (TAUT_TROPOLONE_5/7).
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the walks and ring templates). From each
/// endpoint toward its complementary role: a donor walks
/// single-bond-first (the H leaves along its own single bond), an
/// acceptor walks double/aromatic-bond-first. The 1,5 template
/// terminates on a double bond — a bare single-terminated acceptor
/// is not one (phenol OH's stay fixed) — and needs one endpoint
/// OFF-ring (TAUT_15_NON_RING): two ring-bound endpoints pair only
/// through their own ring templates.
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the TAUT_OTHER block). Four template
/// families, in the reference's order: the 1,5 alt path
/// (TAUT_15_NON_RING), 6-member alt rings (TAUT_4PYRIDINOL_RINGS),
/// pyrazole 5-rings (TAUT_PYRAZOLE_RINGS), and the 1,4 tropolone
/// pattern (TAUT_TROPOLONE_5/7).
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the walks and ring templates). From each
/// endpoint toward its complementary role: a donor walks
/// single-bond-first (the H leaves along its own single bond), an
/// acceptor walks double/aromatic-bond-first. The 1,5 template
/// terminates on a double bond — a bare single-terminated acceptor
/// is not one (phenol OH's stay fixed) — and needs one endpoint
/// OFF-ring (TAUT_15_NON_RING): two ring-bound endpoints pair only
/// through their own ring templates.
/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the TAUT_OTHER block). Four template
/// families, in the reference's order: the 1,5 alt path
/// (TAUT_15_NON_RING), 6-member alt rings (TAUT_4PYRIDINOL_RINGS),
/// pyrazole 5-rings (TAUT_PYRAZOLE_RINGS), and the 1,4 tropolone
/// pattern (TAUT_TROPOLONE_5/7).
fn walk_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut changed = false;

    // TAUT_15_NON_RING: a 1,5 shift along a simple alternating path;
    // at least one of its endpoints off-ring (`nNumAtInRingSystem
    // == 1` — a size-1 ring system IS a non-ring atom).
    for start in 0..n {
        if context.endpoint(start).is_none() {
            continue;
        }
        let mut path = vec![start];
        let mut bonds = Vec::new();
        let mut on_path = vec![false; n];
        on_path[start] = true;
        changed |= walk_alt_path(
            context,
            state,
            &mut path,
            &mut bonds,
            &mut on_path,
            4,
            &mut |state, path, bonds| check_alt_path_15(context, state, path, bonds),
        );
    }

    // TAUT_4PYRIDINOL_RINGS: a ring N (endpoint valence 3, two ring
    // bonds) walks 6-member cycles to the para ring atom bearing an
    // exocyclic =O / -OH / =NH / -NH2; BOTH arcs to the N alternate.
    for start in 0..n {
        if context.rings.size[start] < 6
            || context.neighbors[start].len() != 2
            || !context.is_ring_nitrogen(start)
        {
            continue;
        }
        let mut path = vec![start];
        let mut bonds = Vec::new();
        let mut on_path = vec![false; n];
        on_path[start] = true;
        changed |= walk_ring_path(
            context,
            state,
            &mut path,
            &mut bonds,
            &mut on_path,
            &[],
            5,
            &mut |state, path, bonds| check_6memb_ring(context, state, path, bonds),
        );
    }

    // TAUT_PYRAZOLE_RINGS: two adjacent ring N's (each two ring
    // bonds, one ring system); a fresh pair needs exactly one mobile
    // H between them; the ring walk closes 5-cycles through both.
    for i1 in 0..n {
        if context.rings.size[i1] < 5
            || context.neighbors[i1].len() != 2
            || !context.is_ring_nitrogen(i1)
        {
            continue;
        }
        for &(i2, _) in &context.neighbors[i1] {
            if i2 >= i1 {
                continue;
            }
            if context.rings.system[i2] != context.rings.system[i1] {
                continue;
            }
            if !matches!(
                context.bond_type(state, i1, i2),
                BondType::Single | BondType::Tautom | BondType::Altern
            ) {
                continue;
            }
            if context.neighbors[i2].len() != 2 || !context.is_ring_nitrogen(i2) {
                continue;
            }
            if !state.registered(i1)
                && !state.registered(i2)
                && context.mobile(i1) + context.mobile(i2) != 1
            {
                continue;
            }
            let mut path = vec![i1, i2];
            let mut bonds = vec![context.bond_type(state, i1, i2)];
            let mut on_path = vec![false; n];
            on_path[i1] = true;
            on_path[i2] = true;
            changed |= walk_ring_path(
                context,
                state,
                &mut path,
                &mut bonds,
                &mut on_path,
                &[],
                4,
                &mut |state, path, bonds| check_pyrazole(context, state, path, bonds),
            );
        }
    }

    // TAUT_TROPOLONE_5/7: two ring cut-vertices Q-R (each three
    // neighbors, strict centerpoints, unregistered) joined by a
    // single-like bond, each carrying an EXOCYCLIC endpoint; the H
    // shifts between them across the ring path (the Q-R bond
    // excluded from the alternation).
    for i1 in 0..n {
        if context.rings.size[i1] < 5
            || !context.is_centerpoint_strict(i1)
            || !context.rings.cut_vertex[i1]
            || context.neighbors[i1].len() != 3
            || state.registered(i1)
        {
            continue;
        }
        for &(i2, _) in &context.neighbors[i1] {
            if context.rings.system[i2] != context.rings.system[i1]
                || !context.is_centerpoint_strict(i2)
                || !context.rings.cut_vertex[i2]
                || context.neighbors[i2].len() != 3
                || state.registered(i2)
            {
                continue;
            }
            if !matches!(
                context.bond_type(state, i1, i2),
                BondType::Single | BondType::Tautom | BondType::Altern
            ) {
                continue;
            }
            for &(endpoint1, _) in &context.neighbors[i1] {
                if endpoint1 == i2 {
                    continue;
                }
                let Some(eif1) = context.endpoint(endpoint1) else {
                    continue;
                };
                if context.rings.system[endpoint1] == context.rings.system[i1] {
                    continue;
                }
                let b1 = context.bond_type(state, i1, endpoint1);
                if !matches!(
                    b1,
                    BondType::Single | BondType::Double | BondType::Tautom | BondType::Altern
                ) {
                    continue;
                }
                for &(endpoint2, _) in &context.neighbors[i2] {
                    if endpoint2 == i1 {
                        continue;
                    }
                    let Some(eif2) = context.endpoint(endpoint2) else {
                        continue;
                    };
                    if context.rings.system[endpoint2] == context.rings.system[i2] {
                        continue;
                    }
                    let b2 = context.bond_type(state, i2, endpoint2);
                    if !matches!(
                        b2,
                        BondType::Single | BondType::Double | BondType::Tautom | BondType::Altern
                    ) {
                        continue;
                    }
                    let mut shiftable = false;
                    shiftable |= (state.registered(endpoint1) || context.mobile(endpoint1) > 0)
                        && b1 != BondType::Double
                        && (state.registered(endpoint2) || eif2 == Role::Acceptor)
                        && b2 != BondType::Single;
                    shiftable |= (state.registered(endpoint2) || context.mobile(endpoint2) > 0)
                        && b2 != BondType::Double
                        && (state.registered(endpoint1) || eif1 == Role::Acceptor)
                        && b1 != BondType::Single;
                    if !shiftable {
                        continue;
                    }
                    if b1 == b2 && matches!(b1, BondType::Single | BondType::Double) {
                        continue;
                    }
                    for cycle_bonds in [4, 6] {
                        if cycle_bonds == 6 && context.rings.size[i1] < 7 {
                            continue;
                        }
                        let mut path = vec![i1, i2];
                        let mut bonds = vec![context.bond_type(state, i1, i2)];
                        let mut on_path = vec![false; n];
                        on_path[i1] = true;
                        on_path[i2] = true;
                        changed |= walk_ring_path(
                            context,
                            state,
                            &mut path,
                            &mut bonds,
                            &mut on_path,
                            &[endpoint1, endpoint2],
                            cycle_bonds,
                            &mut |state, path, bonds| {
                                check_tropolone(
                                    context, state, path, bonds, endpoint1, eif1, endpoint2, eif2,
                                )
                            },
                        );
                    }
                }
            }
        }
    }

    changed
}

/// The reference bond classes the templates distinguish. TAUTOM is a
/// bond a previous registration replaced; ALTERN is an aromatic bond.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum BondType {
    Single,
    Double,
    Triple,
    Altern,
    Tautom,
}

impl<'a> Context<'a> {
    fn is_aromatic(&self, a: usize, b: usize) -> bool {
        self.bond_index
            .get(&(a.min(b), a.max(b)))
            .map(|&index| self.aromatic.contains(&index))
            .unwrap_or(false)
    }

    /// The reference's `3 == nGetEndpointInfo(...)`: a classified
    /// nitrogen endpoint.
    fn is_ring_nitrogen(&self, atom: usize) -> bool {
        self.structure.atoms[atom].element == "N" && self.endpoint(atom).is_some()
    }
}

/// `are_alt_bonds`: the path bonds must alternate single/double —
/// ALTERN and TAUTOM match either — and the return is the bond type
/// the final position holds, which the callers read as the pattern
/// bond at the far endpoint. `None` rejects. (The C returns
/// non-single/non-double in the wild-card case, which the callers
/// treat as "no pattern constraint".)
fn are_alt_bonds(bonds: &[BondType]) -> Option<BondType> {
    if bonds.len() < 2 || bonds[0] == BondType::Triple {
        return None;
    }
    let mut expected = match bonds[0] {
        BondType::Single => Some(BondType::Double),
        BondType::Double => Some(BondType::Single),
        _ => None,
    };
    for &bond in &bonds[1..] {
        let wildcard = matches!(bond, BondType::Altern | BondType::Tautom);
        if let Some(want) = expected {
            if bond == want || wildcard {
                expected = Some(if want == BondType::Single {
                    BondType::Double
                } else {
                    BondType::Single
                });
            } else {
                return None;
            }
        } else {
            match bond {
                BondType::Single => expected = Some(BondType::Double),
                BondType::Double => expected = Some(BondType::Single),
                _ if wildcard => {}
                _ => return None,
            }
        }
    }
    Some(match expected {
        Some(BondType::Single) => BondType::Double,
        Some(_) => BondType::Single,
        None => BondType::Altern,
    })
}

/// Mark every marked-able bond of a found path: REPLACE_THE_BOND
/// replaces single, double and alternating bonds with TAUTOM (a
/// triple never; an already-tautom bond stays).
fn mark_path_bonds(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) {
    let ring = bonds.len() == path.len();
    for (step, bond) in path.windows(2).zip(bonds) {
        if matches!(
            *bond,
            BondType::Single | BondType::Double | BondType::Altern
        ) {
            if let Some(&index) = context
                .bond_index
                .get(&(step[0].min(step[1]), step[0].max(step[1])))
            {
                state.marks.insert(index);
            }
        }
    }
    if ring {
        let last = path[path.len() - 1];
        let first = path[0];
        if matches!(
            bonds[bonds.len() - 1],
            BondType::Single | BondType::Double | BondType::Altern
        ) {
            if let Some(&index) = context.bond_index.get(&(first.min(last), first.max(last))) {
                state.marks.insert(index);
            }
        }
    }
}

/// The reference's DFS_FindTautAltPath: simple paths of `target`
/// bonds from the seed; every interior atom a registered endpoint or
/// a strict centerpoint; the far atom is checked by `check` exactly
/// when one endpoint of the shift is off-ring.
fn walk_alt_path(
    context: &Context,
    state: &mut Registration,
    path: &mut Vec<usize>,
    bonds: &mut Vec<BondType>,
    on_path: &mut [bool],
    target: usize,
    check: DfsCheck<'_>,
) -> bool {
    let level = path.len() - 1;
    let atom = path[level];
    let mut changed = false;
    for &(next, _) in &context.neighbors[atom] {
        if level > 0 && next == path[level - 1] {
            continue;
        }
        if on_path[next] {
            continue;
        }
        let bond = context.bond_type(state, atom, next);
        if level + 1 == target {
            if context.rings.size[next] == 1 || context.rings.size[path[0]] == 1 {
                path.push(next);
                bonds.push(bond);
                changed |= check(state, path, bonds);
                bonds.pop();
                path.pop();
            }
        } else if state.registered(next)
            || (context.is_centerpoint_strict(next) && context.has_valence_slack(next))
        {
            on_path[next] = true;
            path.push(next);
            bonds.push(bond);
            changed |= walk_alt_path(context, state, path, bonds, on_path, target, check);
            bonds.pop();
            path.pop();
            on_path[next] = false;
        }
    }
    changed
}

/// The reference's DFS_FindTautInARing: paths around a cycle of
/// `cycle_bonds` bonds closing back into the seed's first atom;
/// interior atoms must be strict centerpoints; forbidden atoms (the
/// 1,4 template's exocyclic endpoints) are never touched.
#[allow(clippy::too_many_arguments)]
fn walk_ring_path(
    context: &Context,
    state: &mut Registration,
    path: &mut Vec<usize>,
    bonds: &mut Vec<BondType>,
    on_path: &mut [bool],
    forbidden: &[usize],
    cycle_bonds: usize,
    check: DfsCheck<'_>,
) -> bool {
    let level = path.len() - 1;
    let atom = path[level];
    let mut changed = false;
    for &(next, _) in &context.neighbors[atom] {
        if forbidden.contains(&next) {
            continue;
        }
        let bond = context.bond_type(state, atom, next);
        if on_path[next] {
            if next == path[0] && level == cycle_bonds {
                bonds.push(bond);
                changed |= check(state, path, bonds);
                bonds.pop();
            }
            continue;
        }
        if !context.is_centerpoint_strict(next) {
            continue;
        }
        if level < cycle_bonds {
            on_path[next] = true;
            path.push(next);
            bonds.push(bond);
            changed |= walk_ring_path(
                context,
                state,
                path,
                bonds,
                on_path,
                forbidden,
                cycle_bonds,
                check,
            );
            bonds.pop();
            path.pop();
            on_path[next] = false;
        }
    }
    changed
}

/// Check15TautPath: both far endpoints heteroatoms; no parity of the
/// path mixes a real single with a real double; each endpoint's role
/// matches its parity; a mobile group must actually shift. The
/// reference's closing `bExistsAnyAltPath` (the BNS cap-aware check)
/// is subsumed here by the walked path itself.
fn check_alt_path_15(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let at1 = path[0];
    let at2 = path[path.len() - 1];
    let Some(eif1) = context.endpoint(at1) else {
        return false;
    };
    let Some(eif2) = context.endpoint(at2) else {
        return false;
    };
    if state.group_of[at1].is_some() && state.group_of[at1] == state.group_of[at2] {
        return false;
    }
    let mut parity = [0u32; 2];
    for (index, bond) in bonds.iter().enumerate() {
        parity[index % 2] |= match bond {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 64,
            BondType::Altern | BondType::Tautom => 4,
        };
    }
    if parity[0] & parity[1] & 3 != 0 || parity[0] & 64 != 0 || parity[1] & 64 != 0 {
        return false;
    }
    if !state.registered(at1) {
        if parity[0] & 1 != 0 && eif1 != Role::Donor {
            return false;
        }
        if parity[0] & 2 != 0 && eif1 != Role::Acceptor {
            return false;
        }
    }
    if !state.registered(at2) {
        if parity[1] & 1 != 0 && eif2 != Role::Donor {
            return false;
        }
        if parity[1] & 2 != 0 && eif2 != Role::Acceptor {
            return false;
        }
    }
    let mut shiftable = false;
    shiftable |= (state.registered(at1) || eif1 == Role::Donor)
        && parity[0] & 2 == 0
        && (state.registered(at2) || eif2 == Role::Acceptor);
    shiftable |= (state.registered(at2) || eif2 == Role::Donor)
        && parity[1] & 2 == 0
        && (state.registered(at1) || eif1 == Role::Acceptor);
    if !shiftable {
        return false;
    }
    state.register(&[at1, at2]);
    mark_path_bonds(context, state, path, bonds);
    true
}

/// Check5MembTautRing (pyrazole): the fresh-pair one-H gate, an
/// alternating 5-cycle, and the closure bond's pattern deciding
/// which N holds the hydrogen.
fn check_pyrazole(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let n1 = path[0];
    let n2 = path[1];
    if context.endpoint(n1).is_none() || context.endpoint(n2).is_none() {
        return false;
    }
    if !state.registered(n1)
        && !state.registered(n2)
        && context.mobile(n1) + context.mobile(n2) != 1
    {
        return false;
    }
    let Some(last) = are_alt_bonds(bonds) else {
        return false;
    };
    match last {
        BondType::Single => {
            if (!state.registered(n1) && context.endpoint(n1) != Some(Role::Donor))
                || (!state.registered(n2) && context.endpoint(n2) != Some(Role::Acceptor))
            {
                return false;
            }
        }
        BondType::Double
            if ((!state.registered(n1) && context.endpoint(n1) != Some(Role::Acceptor))
                || (!state.registered(n2) && context.endpoint(n2) != Some(Role::Donor))) =>
        {
            return false;
        }
        _ => {}
    }
    state.register(&[n1, n2]);
    mark_path_bonds(context, state, path, bonds);
    true
}

/// Check6MembTautRing (4-pyridinol): the para ring atom (three
/// neighbors) carries an exocyclic endpoint; both arcs to the N
/// alternate; a mobile group must actually shift along one.
fn check_6memb_ring(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let n_at = path[0];
    let mid = path[3];
    if context.neighbors[mid].len() != 3 {
        return false;
    }
    let mut o_at = None;
    for &(neighbor, _) in &context.neighbors[mid] {
        if neighbor != path[2] && neighbor != path[4] {
            o_at = Some(neighbor);
            break;
        }
    }
    let Some(o_at) = o_at else {
        return false;
    };
    let exo = context.bond_type(state, mid, o_at);
    if !matches!(
        exo,
        BondType::Single | BondType::Double | BondType::Tautom | BondType::Altern
    ) {
        return false;
    }
    let Some(eif_o) = context.endpoint(o_at) else {
        return false;
    };
    if exo == BondType::Single && !state.registered(o_at) && eif_o != Role::Donor {
        return false;
    }
    let Some(eif_n) = context.endpoint(n_at) else {
        return false;
    };
    let mut shiftable = false;
    shiftable |= (state.registered(o_at) || eif_o == Role::Donor)
        && exo != BondType::Double
        && (state.registered(n_at) || eif_n == Role::Acceptor);
    shiftable |= (state.registered(n_at) || eif_n == Role::Donor)
        && (state.registered(o_at) || eif_o == Role::Acceptor)
        && exo != BondType::Single;
    if !shiftable {
        return false;
    }
    let closure = bonds[bonds.len() - 1];
    if !are_alt_bonds(&[exo, bonds[2], bonds[1], bonds[0]]).is_some()
        || !are_alt_bonds(&[exo, bonds[3], bonds[4], closure]).is_some()
    {
        return false;
    }
    state.register(&[n_at, o_at]);
    mark_path_bonds(context, state, path, bonds);
    true
}

/// Check7MembTautRing (the 1,4 tropolone pattern, its 5-member
/// form): the exocyclic endpoints hang off the seeded ring path;
/// the alternation runs from one exocyclic bond to the other.
#[allow(clippy::too_many_arguments)]
fn check_tropolone(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
    endpoint1: usize,
    eif1: Role,
    endpoint2: usize,
    eif2: Role,
) -> bool {
    let i1 = path[0];
    let i2 = path[1];
    let mut path_bonds: Vec<BondType> = Vec::with_capacity(bonds.len() + 1);
    path_bonds.push(context.bond_type(state, i2, endpoint2));
    path_bonds.extend_from_slice(&bonds[1..]);
    path_bonds.push(context.bond_type(state, i1, endpoint1));
    let Some(last) = are_alt_bonds(&path_bonds) else {
        return false;
    };
    match last {
        BondType::Single => {
            if (!state.registered(endpoint2) && eif2 != Role::Donor)
                || (!state.registered(endpoint1) && eif1 != Role::Acceptor)
            {
                return false;
            }
        }
        BondType::Double
            if ((!state.registered(endpoint2) && eif2 != Role::Acceptor)
                || (!state.registered(endpoint1) && eif1 != Role::Donor)) =>
        {
            return false;
        }
        _ => {}
    }
    state.register(&[endpoint1, endpoint2]);
    let mut marked: Vec<usize> = path.to_vec();
    marked.insert(1, endpoint2);
    marked.push(endpoint1);
    let marked_bonds = vec![BondType::Tautom; marked.len() - 1];
    mark_path_bonds(context, state, &marked, &marked_bonds);
    true
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
            h_pinned: false,
            isotope: 0,
            charge: 0,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
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
    fn kekule_double_bond_anchors_the_pyridinol_pair() {
        // 2-hydroxypyridine skeleton in the kekule form the reference
        // normalizes: the OH oxygen (donor) sits on the anchor carbon
        // whose =N- bond (acceptor) is a real input double, so the
        // 1,3 template pairs them.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms[0] = atom("N");
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 2), // ring N = C1 (the anchor's pi bond)
                (1, 2, 1),
                (2, 3, 2),
                (3, 4, 1),
                (4, 5, 2),
                (5, 0, 1),
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
