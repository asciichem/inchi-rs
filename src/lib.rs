//! Rust InChI engine — the port of the reference engine
//! ([inchi-engine](https://github.com/asciichem/inchi-engine)),
//! gated byte-for-byte by the oracle corpus in
//! [inchi-tests](https://github.com/asciichem/inchi-tests).
//!
//! M1 — the formula layer plus the ±1 proton standardization.

/// A parse or gate failure with a human-readable message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub message: String,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for Error {}

type Result<T> = std::result::Result<T, Error>;

/// Local result alias carrying the borrowed-line lifetime.
type LineResult<'a> = std::result::Result<&'a str, Error>;

fn err<T>(message: impl Into<String>) -> Result<T> {
    Err(Error {
        message: message.into(),
    })
}

const CONTEXT: &str = "inchi::parse_molfile";

#[derive(Debug, Clone, PartialEq)]
pub struct AtomRecord {
    pub element: String,
    /// Explicit H count: MDL code n means n−1 hydrogens (0 = infer).
    pub h_count: u32,
    pub charge: i32,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Structure {
    pub atoms: Vec<AtomRecord>,
    /// V2000 bond types: 1 single, 2 double, 3 triple, 4 aromatic
    /// (counted as 1.5 in valence sums).
    pub bond_types: Vec<(usize, usize, u32)>,
}

/// Parses the V2000 subset the corpus identity molfiles use (the gem
/// writer's output): counts line, atom block, bond block, `M  CHG`
/// property lines.
pub fn parse_molfile(input: &str) -> Result<Structure> {
    let lines: Vec<&str> = input.lines().collect();
    if lines.len() < 5 {
        return err(format!("{CONTEXT}: molfile too short"));
    }
    let counts = lines[3];
    let atom_count: usize = counts
        .get(0..3)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| bad(counts))?;
    let bond_count: usize = counts
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| bad(counts))?;

    let mut atoms = Vec::with_capacity(atom_count);
    for index in 0..atom_count {
        let line = truncated(&lines, 4 + index, "atom block")?;
        let element = line.get(31..34).unwrap_or("").trim().to_string();
        if element.is_empty() {
            return err(format!("{CONTEXT}: atom {index} has no element symbol"));
        }
        let code: u32 = line.get(42..45).unwrap_or("").trim().parse().unwrap_or(0);
        atoms.push(AtomRecord {
            element,
            h_count: code.saturating_sub(1),
            charge: 0,
        });
    }

    let mut bond_types = Vec::with_capacity(bond_count);
    for index in 0..bond_count {
        let line = truncated(&lines, 4 + atom_count + index, "bond block")?;
        let from: usize = line
            .get(0..3)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        let to: usize = line
            .get(3..6)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        let kind: u32 = line
            .get(6..9)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        bond_types.push((from - 1, to - 1, kind));
    }

    for line in &lines[4 + atom_count + bond_count..] {
        let is_charge = line.starts_with("M  CHG");
        if !(is_charge || line.starts_with("M  ISO")) {
            continue;
        }
        let values: Vec<&str> = line[6..].split_whitespace().collect();
        for pair in values[1..].chunks(2) {
            if pair.len() != 2 {
                continue;
            }
            let index: usize = pair[0].parse().map_err(|_| bad(line))?;
            let value: i32 = pair[1].parse().map_err(|_| bad(line))?;
            let atom = atoms.get_mut(index - 1).ok_or_else(|| bad(line))?;
            if is_charge {
                atom.charge = value;
            }
        }
    }

    Ok(Structure { atoms, bond_types })
}

fn truncated<'a>(lines: &[&'a str], index: usize, block: &str) -> LineResult<'a> {
    lines.get(index).copied().ok_or_else(|| Error {
        message: format!("{CONTEXT}: truncated {block}"),
    })
}

fn bad(line: &str) -> Error {
    Error {
        message: format!("{CONTEXT}: malformed line: {line:?}"),
    }
}

fn bond_order(kind: u32) -> f64 {
    match kind {
        2 => 2.0,
        3 => 3.0,
        4 => 1.5,
        _ => 1.0,
    }
}

// Standard valences of the organic subset; hydrogens implied by the
// gap between these and the bond-order sum.
fn valence(element: &str) -> f64 {
    match element {
        "C" => 4.0,
        "N" => 3.0,
        "O" => 2.0,
        "S" => 6.0,
        "P" => 5.0,
        "F" | "Cl" | "Br" | "I" | "H" => 1.0,
        _ => 0.0,
    }
}

impl Structure {
    pub fn total_charge(&self) -> i32 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// Element → count, with implicit hydrogens: MDL H-count codes
    /// win, otherwise the valence gap (the oracle does the same for
    /// code-0 atoms — a bare `O` line derives H2O). Anionic
    /// heteroatoms take zero: their standardization crosses the
    /// `/p` layer.
    pub fn composition(&self) -> Result<Vec<(String, u32)>> {
        use std::collections::BTreeMap;
        let mut sums = vec![0.0; self.atoms.len()];
        for &(from, to, kind) in &self.bond_types {
            let order = bond_order(kind);
            let first = sums.get_mut(from).ok_or_else(|| bad("bond index"))?;
            *first += order;
            let second = sums.get_mut(to).ok_or_else(|| bad("bond index"))?;
            *second += order;
        }
        let mut map: BTreeMap<String, u32> = BTreeMap::new();
        for (atom, sum) in self.atoms.iter().zip(&sums) {
            *map.entry(atom.element.clone()).or_insert(0) += 1;
            let hydrogens = if atom.h_count > 0 {
                atom.h_count
            } else if atom.charge < 0 {
                0
            } else {
                (valence(&atom.element) - sum).max(0.0) as u32
            };
            if hydrogens > 0 {
                *map.entry("H".to_string()).or_insert(0) += hydrogens;
            }
        }
        Ok(map.into_iter().collect())
    }

    /// The `/formula` layer plus the ±1 proton standardization the
    /// charged corpus cases carry: an anion is recorded as its
    /// neutral parent with one MORE hydrogen (`/p-1`), a cation with
    /// one LESS (`/p+1`). Proton layers are signed strings.
    pub fn formula_layer(&self) -> Result<(String, Option<String>)> {
        let mut composition = self.composition()?;
        let charge = self.total_charge();
        let proton = match charge {
            0 => None,
            -1 => {
                add_hydrogen(&mut composition, 1);
                Some("-1".to_string())
            }
            1 => {
                let has_h = composition.iter().any(|(e, c)| e == "H" && *c > 0);
                if !has_h {
                    return err("cation without hydrogen to standardize");
                }
                add_hydrogen(&mut composition, -1);
                Some("+1".to_string())
            }
            other => {
                return err(format!(
                    "{CONTEXT}: |charge| = {other} beyond the ±1 corpus gate"
                ))
            }
        };
        Ok((hill_string(&composition), proton))
    }
}

fn add_hydrogen(composition: &mut Vec<(String, u32)>, delta: i32) {
    for entry in composition.iter_mut() {
        if entry.0 == "H" {
            entry.1 = (entry.1 as i64 + delta as i64).max(0) as u32;
            return;
        }
    }
    if delta > 0 {
        composition.push(("H".to_string(), delta as u32));
    }
}

// Hill order with the InChI twist: carbon first (when present), then
// hydrogen, then the rest alphabetically; counts omitted when 1.
// Without carbon: pure alphabetical.
fn hill_string(composition: &[(String, u32)]) -> String {
    let has_carbon = composition.iter().any(|(e, _)| e == "C");
    let rank = |element: &str| match (has_carbon, element) {
        (true, "C") => 0,
        (true, "H") => 1,
        (true, _) => 2,
        (false, _) => 0,
    };
    let mut entries: Vec<&(String, u32)> = composition.iter().collect();
    entries.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)).then(a.0.cmp(&b.0)));
    entries
        .iter()
        .filter(|(_, count)| *count > 0)
        .map(|(element, count)| {
            if *count > 1 {
                format!("{element}{count}")
            } else {
                element.clone()
            }
        })
        .collect()
}
