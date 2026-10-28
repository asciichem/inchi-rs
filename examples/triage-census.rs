//! Triage census: for every triaged fixture, compare the assembled
//! identifier with the oracle recording — the retirement radar.
use std::fs;

use inchi_rs::Engine;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    use std::path::PathBuf;
    let triage: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("corpus/triage.json")).expect("triage"))
            .expect("valid");
    let mut files: Vec<PathBuf> = fs::read_dir(dir.join("corpus/fixtures"))
        .expect("dir")
        .filter_map(|e| {
            let p = e.expect("entry").path();
            let n = p.file_name()?.to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let mut passing = Vec::new();
    let mut failing_triaged = Vec::new();
    let mut total = 0usize;
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("json");
        for case in parsed {
            let id = case["id"].as_str().expect("id").to_string();
            if triage.get(&id).is_none() {
                continue;
            }
            total += 1;
            let Ok(engine) = Engine::from_molfile(case["molfile"].as_str().unwrap()) else {
                continue;
            };
            let Ok(assembled) = engine.inchi_string() else {
                continue;
            };
            if assembled == case["inchi"].as_str().unwrap_or("") {
                passing.push(id);
            } else {
                failing_triaged.push((
                    id,
                    assembled,
                    case["inchi"].as_str().unwrap_or("").to_string(),
                ));
            }
        }
    }
    println!("triaged {total}, now passing {}:", passing.len());
    for id in &passing {
        println!("  {id}");
    }
    println!("triaged, still failing {}:", failing_triaged.len());
    for (id, assembled, oracle) in &failing_triaged {
        println!("  {id}\n    ours   {assembled}\n    oracle {oracle}");
    }
}
