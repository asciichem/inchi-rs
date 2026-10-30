//! The engine façade: parses a molfile and runs the registered
//! layers. Closed for modification, open for extension — new layers
//! join the registry in `layers::default_layers` without this file
//! changing.

use crate::error::{err, Result};
use crate::graph::Structure;
use crate::layers::{default_layers, standardization, LayerEmitter, MANDATORY_KINDS};
use crate::molfile;

pub struct Engine {
    structure: Structure,
    layers: Vec<Box<dyn LayerEmitter>>,
}

impl Engine {
    /// Parses the V2000 molfile subset, normalizes internal salts,
    /// folds explicit hydrogen atoms onto their carriers, and equips
    /// the default layer registry.
    pub fn from_molfile(input: &str) -> Result<Self> {
        let mut structure = molfile::parse(input)?;
        standardization::normalize_internal_salts(&mut structure);
        structure.fold_explicit_hydrogens();
        // The (DE)PROTONATION normalization: charged structures read
        // as their conjugate form before any layer runs; the pinned
        // H/charge state is final.
        let mut proton_shift = 0;
        if let Some(std) = crate::canonical::protons::standardize(&structure) {
            proton_shift = std.removed_protons as i32 - std.added_protons as i32;
            for (atom, h) in structure.atoms.iter_mut().zip(&std.h) {
                atom.h_count = *h;
                atom.h_pinned = true;
            }
            for (atom, charge) in structure.atoms.iter_mut().zip(&std.charge) {
                atom.charge = *charge;
            }
        }
        proton_shift += remove_group_protons(&mut structure);
        Ok(Engine {
            structure,
            layers: default_layers(proton_shift),
        })
    }

    /// The normalized structure — diagnostics for the triage tools.
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// The derived segments, tagged by layer kind — everything the
    /// registered layers emit. Milestones add layers; the gates match
    /// tagged segments against the oracle's recording by kind.
    pub fn segments(&self) -> Result<Vec<(Option<char>, String)>> {
        let mut out = Vec::new();
        for layer in &self.layers {
            if let Some(segment) = layer.emit(&self.structure)? {
                out.push((layer.segment_kind(), segment));
            }
        }
        Ok(out)
    }

    /// Assembles the full InChI string, segments sorted by declared
    /// layer order. Errors until the registered layers cover the
    /// mandatory kinds (formula + connection): a partial layer set
    /// must never present itself as an identifier.
    pub fn inchi_string(&self) -> Result<String> {
        let registered: Vec<Option<char>> = self.layers.iter().map(|l| l.segment_kind()).collect();
        for mandatory in &MANDATORY_KINDS {
            if !registered.contains(mandatory) {
                return err(format!(
                    "inchi_string: mandatory layer kind {:?} not registered yet — \
                     the engine's layer set is milestone-partial",
                    mandatory
                        .map(|c| c.to_string())
                        .unwrap_or_else(|| "formula".to_string())
                ));
            }
        }

        let mut sorted: Vec<&Box<dyn LayerEmitter>> = self.layers.iter().collect();
        sorted.sort_by_key(|layer| layer.order());
        let mut out = String::from("InChI=1S/");
        let mut appended = false;
        for layer in sorted {
            if let Some(segment) = layer.emit(&self.structure)? {
                if appended {
                    out.push('/');
                }
                out.push_str(&segment);
                appended = true;
            }
        }
        Ok(out)
    }

    /// The InChIKey over the assembled identifier — the truncated
    /// SHA-256 key, computed natively (see [`crate::key`]).
    pub fn inchikey(&self) -> Result<String> {
        crate::key::inchikey(&self.inchi_string()?)
    }
}

/// The hard group-H removal (the reference's HardRemoveAcidicProtons,
/// TECHMAN 5.2b): a cationic structure whose registered mobile-H
/// group still holds a removable proton sheds it onto the group's
/// acidic oxygen as a fixed −1 — the cation stays put, the pinned
/// total reads neutral, and the shift surfaces as `/p+1` with no `/q`
/// and no mobile group. Each removal is gated by the reference's
/// REM_PROTON route test over the two virtual group stars (the
/// H-donor star and the acidic =O/N-O star): only a proton whose
/// donor-to-acceptor skeleton route the balanced search can cover
/// leaves the structure.
fn remove_group_protons(structure: &mut Structure) -> i32 {
    let mut removals = 0;
    loop {
        // The quota: one removal per unit of positive charge. The +2
        // family lands when the charged super-group merge follows.
        let total = structure.total_charge() + removals;
        if total <= 0 {
            break;
        }
        // The fixed cation the shift leaves behind: an N/P(+) with no
        // hydrogen of its own (the quaternary centers included).
        let plus = structure.atoms.iter().position(|atom| {
            atom.charge == 1
                && atom.h_count == 0
                && matches!(
                    crate::canonical::protons::el_number(&atom.element),
                    Some(7 | 15)
                )
        });
        let Some(plus) = plus else {
            break;
        };
        // The derived hydrogen counts (implicit H never shows in
        // `h_count`), the canonical ranks and the registered mobile
        // groups over them.
        let hydrogens = structure.hydrogens();
        let ranks = crate::layers::hydrogen::canonical_ranks_for(structure);
        let groups = crate::canonical::tautomer::mobile_groups(structure, &hydrogens, &ranks);
        let mut membership: Vec<Option<usize>> = vec![None; structure.atoms.len()];
        for (index, group) in groups.iter().enumerate() {
            for &rank in &group.ranks {
                if let Some(atom) = ranks.iter().position(|&r| r == rank) {
                    membership[atom] = Some(index);
                }
            }
        }
        // The holder: the group member still carrying its hydrogen.
        let Some(holder) = (0..structure.atoms.len())
            .find(|&atom| membership[atom].is_some() && hydrogens[atom] > 0)
        else {
            break;
        };
        // The acceptor: the group's double-bonded, hydrogen-free
        // member — where the moved proton turns into a fixed charge.
        let sums = structure.bond_order_sums();
        let degrees = structure.degrees();
        let members: Vec<usize> = groups[membership[holder].unwrap()]
            .ranks
            .iter()
            .filter_map(|&rank| ranks.iter().position(|&r| r == rank))
            .collect();
        let acceptor = members
            .iter()
            .copied()
            .find(|&atom| hydrogens[atom] == 0 && sums[atom] > degrees[atom] as f64)
            .unwrap_or(holder);
        let _ = plus;
        // The route: donor star vs acidic star over the whole
        // structure — no route, no removal.
        let route = rem_proton_route_for(structure, &hydrogens, holder, acceptor);
        if std::env::var_os("INCHI_TRACE_PROTON").is_some() {
            eprintln!("GROUP-REM plus={plus} holder={holder} acceptor={acceptor} route={route}");
        }
        if !route {
            break;
        }
        // Apply: freeze the derived counts (the ±1 heuristics must not
        // re-read the half-moved state), then move the proton.
        if !structure.atoms.iter().all(|atom| atom.h_pinned) {
            for (atom, h) in structure.atoms.iter_mut().zip(&hydrogens) {
                atom.h_count = *h;
                atom.h_pinned = true;
            }
        }
        structure.atoms[holder].h_count = hydrogens[holder] - 1;
        structure.atoms[acceptor].charge -= 1;
        removals += 1;
    }
    removals
}

/// The REM_PROTON route test over a fresh network for this pass: one
/// t-group star over the H-donor heteroatoms, one over the acidic
/// =O/N-O oxygens, both as the reference builds them for the hard
/// removal (AR_HARD classes, neutral-valence gate, member-edge
/// capacity raises).
fn rem_proton_route_for(
    structure: &Structure,
    hydrogens: &[u32],
    holder: usize,
    acceptor: usize,
) -> bool {
    let context = crate::canonical::protons::ProtonContext::new(structure);
    let (donors, acids) = crate::canonical::protons::hard_removal_stars(&context);
    let n = structure.atoms.len();
    let mut bond_adjacency: Vec<Vec<(usize, i32, bool)>> = vec![Vec::new(); n];
    let mut single_like: Vec<Vec<bool>> = vec![Vec::new(); n];
    for &(from, to, kind) in &structure.bonds {
        let flow = match kind {
            2 => 1,
            3 => 2,
            _ => 0,
        };
        bond_adjacency[from].push((to, flow, false));
        bond_adjacency[to].push((from, flow, false));
        single_like[from].push(true);
        single_like[to].push(true);
    }
    let max_flow: Vec<i32> = (0..n)
        .map(|atom| (context.cbv2[atom] / 2 - context.neighbors[atom].len() as i32).max(0))
        .collect();
    let group_of: Vec<Option<usize>> = (0..n)
        .map(|atom| {
            if donors[atom] {
                Some(0)
            } else if acids[atom] {
                Some(1)
            } else {
                None
            }
        })
        .collect();
    // A star member pushes only when it actually holds hydrogen (the
    // reference's member-edge flow is the capped H count).
    let donor: Vec<bool> = (0..n)
        .map(|atom| donors[atom] && hydrogens[atom] > 0)
        .collect();
    let net = crate::canonical::bns::Network::build(&crate::canonical::bns::Build {
        bond_adjacency: &bond_adjacency,
        max_flow: &max_flow,
        group_of: &group_of,
        donor: &donor,
        single_like: &single_like,
        c_point_of: &[None; 0],
        c_charged: &[false; 0],
    });
    crate::canonical::bns::rem_proton_route(&net, holder, acceptor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Structure;
    use crate::layers::order;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    struct FixedSegment {
        order: u8,
        kind: Option<char>,
        text: &'static str,
    }

    impl LayerEmitter for FixedSegment {
        fn name(&self) -> &'static str {
            "fixed"
        }

        fn order(&self) -> u8 {
            self.order
        }

        fn segment_kind(&self) -> Option<char> {
            self.kind
        }

        fn emit(&self, _structure: &Structure) -> Result<Option<String>> {
            Ok(Some(self.text.to_string()))
        }
    }

    #[test]
    fn a_milestone_partial_layer_set_refuses_to_assemble_an_inchi() {
        // No connection layer registered: the set is milestone-partial
        // and must not present a prefix as an InChI.
        let engine = Engine {
            structure: crate::parse_molfile(WATER).expect("parses"),
            layers: vec![Box::new(FixedSegment {
                order: order::FORMULA,
                kind: None,
                text: "H2O",
            })],
        };
        let error = engine.inchi_string().unwrap_err();
        assert!(error.message().contains("mandatory layer kind"));
    }

    #[test]
    fn water_assembles_from_the_registered_layers() {
        let engine = Engine::from_molfile(WATER).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn assembly_follows_declared_layer_order_not_registry_order() {
        // Registry scrambled on purpose; the assembler must sort by
        // `order()` and join only the segments that emit.
        let engine = Engine {
            structure: Structure::default(),
            layers: vec![
                Box::new(FixedSegment {
                    order: order::P,
                    kind: Some('p'),
                    text: "p+1",
                }),
                Box::new(FixedSegment {
                    order: order::CONNECTION,
                    kind: Some('c'),
                    text: "c1-2",
                }),
                Box::new(FixedSegment {
                    order: order::FORMULA,
                    kind: None,
                    text: "H2O",
                }),
            ],
        };
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/c1-2/p+1");
    }
}
