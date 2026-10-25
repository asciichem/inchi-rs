//! Triage diagnostics: dump our emitted segments for fixture ids
//! given as args (substring match), with the recorded oracle line
//! alongside — the compare-always corpus gate's companion.
use std::fs;
use std::path::PathBuf;

use inchi_rs::Engine;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let mut files: Vec<PathBuf> = fs::read_dir(dir.join("corpus/fixtures"))
        .expect("dir")
        .filter_map(|e| {
            let p = e.expect("entry").path();
            let n = p.file_name()?.to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let want: Vec<String> = std::env::args().skip(1).collect();
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("json");
        for case in parsed {
            let id = case["id"].as_str().expect("id").to_string();
            if want.iter().any(|w| id.contains(w)) {
                let engine = Engine::from_molfile(case["molfile"].as_str().unwrap()).unwrap();
                println!("{id}:");
                if let Ok(segs) = engine.segments() {
                    for (kind, seg) in segs {
                        println!("  {:?} {seg}", kind);
                    }
                }
                println!("  oracle: {}", case["inchi"].as_str().unwrap());
            }
        }
    }
}
