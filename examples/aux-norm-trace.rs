//! Trace exactly what aux.rs sees for the norm digit: the working
//! structure's groups, charge plan, pinned state.
use std::fs;

use inchi_rs::Engine;

fn main() {
    let id = std::env::args().nth(1).expect("id");
    for file in fs::read_dir("../inchi-tests/corpus/fixtures").expect("dir") {
        let p = file.expect("e").path();
        let n = p.file_name().unwrap().to_string_lossy().to_string();
        if !n.starts_with("identity") || !n.ends_with(".json") {
            continue;
        }
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&p).expect("read")).expect("json");
        for case in cases {
            if case["id"].as_str().unwrap() == format!("identity/{id}") {
                // Rebuild the PRE-standardization structure the way
                // Engine::from_molfile does before protons::standardize.
                let mut pre = inchi_rs::parse_molfile(case["molfile"].as_str().unwrap()).unwrap();
                inchi_rs::layers::standardization::normalize_internal_salts(&mut pre);
                pre.fold_explicit_hydrogens();
                let pre_counts = pre.hydrogens();
                let pre_ranks = inchi_rs::layers::hydrogen::canonical_ranks_for(&pre);
                let pre_groups =
                    inchi_rs::canonical::tautomer::mobile_groups(&pre, &pre_counts, &pre_ranks);
                println!("PRE groups:");
                for g in &pre_groups {
                    let members: Vec<String> = g.ranks.iter().map(|r| r.to_string()).collect();
                    println!("  group count={} charge={}: ranks {}", g.count, g.charge, members.join(","));
                }
                let engine = Engine::from_molfile(case["molfile"].as_str().unwrap()).unwrap();
                let s = engine.structure();
                let counts = inchi_rs::layers::standardization::adjusted_hydrogens(s);
                let ranks = inchi_rs::layers::hydrogen::canonical_ranks_for(s);
                let groups = inchi_rs::canonical::tautomer::mobile_groups(s, &counts, &ranks);
                let pinned = s.atoms.iter().all(|a| a.h_pinned);
                println!("pinned={pinned} protons_hard_moved={} total_charge={}", s.protons_hard_moved, s.total_charge());
                for (i, a) in s.atoms.iter().enumerate() {
                    if a.charge != 0 || counts[i] > 0 && matches!(a.element.as_str(), "O"|"N"|"S"|"Se"|"Te") {
                        println!(
                            "  atom {i}: {} charge={:+} h={} pinned={} rank={} counts={}",
                            a.element, a.charge, a.h_count, a.h_pinned, ranks[i], counts[i]
                        );
                    }
                }
                for g in &groups {
                    let members: Vec<String> = g.ranks.iter().map(|r| r.to_string()).collect();
                    println!("group count={} charge={}: ranks {}", g.count, g.charge, members.join(","));
                }
                let plan = inchi_rs::layers::standardization::charge_plan(s);
                println!("charge_plan: {:?}", plan);
                println!("inchi: {}", engine.inchi_string().unwrap());
                println!("aux:   {}", engine.aux_string().unwrap());
                return;
            }
        }
    }
}
