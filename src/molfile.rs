//! The V2000 molfile reader — the identity pipeline's input edge.
//! Encapsulated on purpose: callers see a `Structure`, never the
//! column arithmetic; only the subset the corpus identity molfiles
//! use (the gem writer's output) is supported, and unsupported
//! constructs are reported rather than guessed.

use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

const CONTEXT: &str = "molfile";

/// Reads the V2000 subset: counts line, atom block, bond block, and
/// `M  CHG` / `M  ISO` property lines. `M  ISO` is consumed and
/// ignored at this level — isotopes do not affect any M1 layer.
pub fn parse(input: &str) -> Result<Structure> {
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
        .map_err(|_| bad_counts(counts))?;
    let bond_count: usize = counts
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| bad_counts(counts))?;

    let mut atoms = Vec::with_capacity(atom_count);
    for index in 0..atom_count {
        let line = line_or(&lines, 4 + index, "atom block")?;
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

    let mut bonds = Vec::with_capacity(bond_count);
    for index in 0..bond_count {
        let line = line_or(&lines, 4 + atom_count + index, "bond block")?;
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
        bonds.push((from - 1, to - 1, kind));
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

    Ok(Structure { atoms, bonds })
}

fn line_or<'a>(lines: &'a [&'a str], index: usize, block: &str) -> Result<&'a str> {
    lines
        .get(index)
        .copied()
        .ok_or_else(|| Error::new(format!("{CONTEXT}: truncated {block}")))
}

use crate::error::Error;

fn bad_counts(line: &str) -> Error {
    Error::new(format!("{CONTEXT}: malformed counts line: {line:?}"))
}

fn bad(line: &str) -> Error {
    Error::new(format!("{CONTEXT}: malformed line: {line:?}"))
}
