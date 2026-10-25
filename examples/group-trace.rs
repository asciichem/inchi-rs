//! Trace mobile groups: members with elements.
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
                let engine = Engine::from_molfile(case["molfile"].as_str().unwrap()).unwrap();
                let s = engine.structure();
                let counts = inchi_rs::layers::standardization::adjusted_hydrogens(s);
                let groups = {
                    let (adjacency, invariants) = inchi_rs::canonical::graph_of(s);
                    let mask = inchi_rs::canonical::tautomer::endpoint_mask(s, &counts);
                    let ranks = inchi_rs::canonical::canonical_ranks(
                        &adjacency,
                        &invariants,
                        &inchi_rs::canonical::tautomer::without_mobile_h(&counts, &mask),
                    );
                    inchi_rs::canonical::tautomer::mobile_groups(s, &counts, &ranks)
                };
                let mask = inchi_rs::canonical::tautomer::endpoint_mask(s, &counts);
                for (i, a) in s.atoms.iter().enumerate() {
                    if mask[i] {
                        println!("masked atom {i}: {} h={}", a.element, counts[i]);
                    }
                }
                for g in &groups {
                    let members: Vec<String> = g.ranks.iter().map(|r| format!("{}", r)).collect();
                    println!("group H{}: ranks {}", g.count, members.join(","));
                }
                println!(
                    "atoms: {}",
                    s.atoms
                        .iter()
                        .map(|a| a.element.as_str())
                        .collect::<Vec<_>>()
                        .join(",")
                );
                return;
            }
        }
    }
}
