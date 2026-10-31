use chemomer::Engine;
use std::fs;

fn main() {
    let dir = "../inchi-tests/corpus/fixtures";
    let mut files: Vec<_> = fs::read_dir(dir)
        .expect("dir")
        .filter_map(|e| {
            let p = e.expect("e").path();
            let n = p.file_name().unwrap().to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let mut mism = 0;
    let mut total = 0;
    for file in files {
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("json");
        for case in cases {
            let Some(aux) = case["aux"].as_str() else {
                continue;
            };
            let oracle_norm: i32 = aux
                .split("AuxInfo=1/")
                .nth(1)
                .and_then(|rest| rest.split('/').next())
                .and_then(|d| d.parse().ok())
                .unwrap_or(-1);
            let engine = match Engine::from_molfile(case["molfile"].as_str().unwrap()) {
                Ok(e) => e,
                Err(_) => continue,
            };
            let s = engine.structure();
            let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
            let ranks = chemomer::layers::hydrogen::canonical_ranks_for(s);
            let groups = chemomer::canonical::tautomer::mobile_groups(s, &counts, &ranks);
            let real_h = |g: &chemomer::canonical::tautomer::MobileGroup| g.count > g.charge;
            let taut = groups.iter().any(|g| real_h(g));
            let minus = groups.iter().any(|g| real_h(g) && g.charge > 0);
            let hard = s.protons_hard_moved;
            let pinned = s.atoms.iter().all(|a| a.h_pinned);
            let has_p = if pinned {
                engine.proton_shift() != 0
            } else {
                chemomer::layers::standardization::charge_plan(s).p != 0
            };
            let iso_mob = exchangeable(s, &counts);
            let bit0 = (taut || has_p || iso_mob) as i32;
            let salt = ((minus && !hard) || (!minus && hard && taut)) as i32;
            let pred = if bit0 == 0 {
                0
            } else {
                bit0 + salt + 4 * hard as i32
            };
            total += 1;
            if pred != oracle_norm {
                mism += 1;
                if mism <= 10 {
                    println!("MISM {} oracle={oracle_norm} pred={pred} taut={taut} minus={minus} hard={hard} p={has_p} iso={iso_mob}", case["id"].as_str().unwrap());
                }
            }
        }
    }
    println!("census: {mism} mismatches of {total}");
}

fn exchangeable(s: &chemomer::graph::Structure, counts: &[u32]) -> bool {
    let isotopic = !s.isotopic_h.is_empty() || s.atoms.iter().any(|a| a.isotope != 0);
    isotopic
        && s.atoms
            .iter()
            .enumerate()
            .any(|(i, a)| a.element != "C" && a.element != "H" && counts[i] > 0)
}
