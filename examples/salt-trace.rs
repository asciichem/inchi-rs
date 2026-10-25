//! Trace the annihilation for one fixture.
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
                for (i, a) in engine.structure().atoms.iter().enumerate() {
                    if a.charge != 0 {
                        println!("atom {i} {} charge {}", a.element, a.charge);
                    }
                }
                println!(
                    "plan p: {:?}",
                    inchi_rs::layers::standardization::proton_layer(engine.structure())
                );
                println!(
                    "plan q: {:?}",
                    inchi_rs::layers::standardization::charge_layer(engine.structure())
                );
                return;
            }
        }
    }
    println!("not found");
}
