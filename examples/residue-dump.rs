//! M3 diagnostics: for every triaged fixture, our /h groups vs the
//! oracle's, with the delta classified.

use std::fs;
use std::path::PathBuf;

use inchi_rs::Engine;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let fixtures = dir.join("corpus/fixtures");
    let triage: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("corpus/triage.json")).expect("triage"))
            .expect("valid triage");
    let mut files: Vec<PathBuf> = fs::read_dir(&fixtures)
        .expect("dir")
        .filter_map(|e| {
            let p = e.expect("entry").path();
            let n = p.file_name()?.to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("json");
        for case in parsed {
            let id = case["id"].as_str().expect("id").to_string();
            if triage.get(&id).is_none() {
                continue;
            }
            let engine = Engine::from_molfile(case["molfile"].as_str().unwrap()).unwrap();
            let segments = engine.segments().unwrap();
            let h = segments
                .iter()
                .find(|(kind, _)| *kind == Some('h'))
                .map(|(_, s)| s.clone())
                .unwrap_or_default();
            let recorded = case["inchi"].as_str().unwrap();
            let oracle_h = recorded
                .split('/')
                .find(|s| s.starts_with('h'))
                .map(|s| s[1..].to_string())
                .unwrap_or_default();
            if h == oracle_h {
                continue; // /h matches; only other layers triaged
            }
            let groups = |s: &str| -> Vec<(String, Vec<String>)> {
                let mut out = Vec::new();
                let mut chars = s.chars().peekable();
                while let Some(c) = chars.next() {
                    if c == '(' {
                        let mut inner = String::new();
                        for d in chars.by_ref() {
                            if d == ')' {
                                break;
                            }
                            inner.push(d);
                        }
                        let mut parts = inner.split(',');
                        let head = parts.next().unwrap_or("").to_string();
                        let ranks: Vec<String> = parts.map(|r| r.to_string()).collect();
                        out.push((head, ranks));
                    }
                }
                out
            };
            let og = groups(&h);
            let rg = groups(&oracle_h);
            let key = |head: &str| head.trim_start_matches('H').to_string();
            let osets: Vec<(String, std::collections::BTreeSet<String>)> = og
                .iter()
                .map(|(head, ranks)| (key(head), ranks.iter().cloned().collect()))
                .collect();
            let rsets: Vec<(String, std::collections::BTreeSet<String>)> = rg
                .iter()
                .map(|(head, ranks)| (key(head), ranks.iter().cloned().collect()))
                .collect();
            let kind = if osets.len() > rsets.len() {
                "extra-groups"
            } else if rsets.len() > osets.len() {
                "missing-groups"
            } else {
                "members-differ"
            };
            println!("{id}\tours h{h}\toracle h{oracle_h}\t{kind}");
        }
    }
}
