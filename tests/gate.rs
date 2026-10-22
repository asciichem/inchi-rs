//! The M1 gate: the same inchi-tests corpus the Ruby reference
//! passes — the derived formula layer must equal the recorded oracle
//! segment for every identity fixture, 100%, no skips.

use std::fs;
use std::path::PathBuf;

fn corpus() -> Vec<serde_json::Value> {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let path = dir.join("corpus/fixtures/identity.json");
    serde_json::from_str(&fs::read_to_string(&path).expect("corpus corpus/fixtures/identity.json"))
        .expect("valid fixtures")
}

#[test]
fn formula_layer_agrees_with_the_oracle() {
    let cases = corpus();
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let molfile = case["molfile"].as_str().expect("molfile");
        let recorded = case["inchi"].as_str().expect("inchi");
        let expected = recorded.split('/').nth(1).expect("formula segment");

        match inchi_rs::parse_molfile(molfile).and_then(|s| s.formula_layer()) {
            Ok((formula, proton)) => {
                if formula != expected {
                    failures.push(format!("{id}: formula {formula:?} != oracle {expected:?}"));
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
            }
            Err(error) => failures.push(format!("{id}: {error}")),
        }
    }

    assert!(
        failures.is_empty(),
        "{} formula-layer failure(s):\n{}\n(pass: {}/{} cases)",
        failures.len(),
        failures.join("\n"),
        cases.len() - failures.len(),
        cases.len()
    );
}
