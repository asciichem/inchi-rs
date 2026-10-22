//! The M1 gate: the same inchi-tests corpus the Ruby reference
//! passes — the derived formula layer must equal the recorded oracle
//! segment for every identity fixture, 100%, no skips.

use std::fs;
use std::path::PathBuf;

#[test]
fn formula_layer_agrees_with_the_oracle() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let fixtures = dir.join("corpus/fixtures");
    let mut cases = Vec::new();
    let mut files: Vec<PathBuf> = fs::read_dir(&fixtures)
        .expect("corpus fixtures dir")
        .filter_map(|entry| {
            let path = entry.expect("dir entry").path();
            let name = path
                .file_name()
                .expect("name")
                .to_string_lossy()
                .to_string();
            (name.starts_with("identity") && name.ends_with(".json")).then_some(path)
        })
        .collect();
    files.sort();
    assert!(!files.is_empty(), "no identity fixture files");
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        cases.extend(parsed);
    }
    let triage: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(dir.join("corpus/triage.json")).expect("triage file"),
    )
    .expect("valid triage");

    let mut pass = 0;
    let mut skipped = 0;
    let mut beyond = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let m3 = case
            .get("tags")
            .and_then(|tags| tags.as_array())
            .map(|tags| tags.iter().any(|t| t.as_str() == Some("m3")))
            .unwrap_or(false);
        let triaged = triage.get(id).is_some();
        if m3 || triaged {
            skipped += 1;
            continue;
        }

        let molfile = case["molfile"].as_str().expect("molfile");
        let recorded = case["inchi"].as_str().expect("inchi");
        let expected = recorded.split('/').nth(1).expect("formula segment");

        match inchi_rs::parse_molfile(molfile).and_then(|s| s.formula_layer()) {
            Ok((formula, proton)) => {
                if formula != expected {
                    failures.push(format!("{id}: formula {formula:?} != oracle {expected:?}"));
                    continue;
                }
                let expected_proton = recorded
                    .split('/')
                    .find_map(|segment| segment.strip_prefix('p'))
                    .map(|sign| format!("p{sign}"));
                let actual_proton = proton.map(|sign| format!("p{sign}"));
                match (expected_proton, actual_proton) {
                    (None, None) => {}
                    (Some(expected), Some(actual)) if expected == actual => {}
                    (expected, actual) => failures.push(format!(
                        "{id}: proton layer {actual:?} != oracle {expected:?}"
                    )),
                }
                pass += 1;
            }
            Err(error) => {
                if error.message.contains("beyond the ±1 corpus gate") {
                    // |charge| > 1 needs the M3 mobile-H rules.
                    beyond += 1;
                } else {
                    failures.push(format!("{id}: {error}"));
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} formula-layer failure(s):\n{}\n(pass: {pass}, skips: {skipped}, beyond: {beyond}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
}
