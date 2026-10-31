//! The gate: the same inchi-tests corpus the reference engine passes.
//! Contract: FULL byte equality — the engine's assembled identifier
//! equals the oracle's recording exactly, every layer present in the
//! oracle's order (formula, /c, /h, /q, /p, stereo, /i), plus the
//! native InChIKey. Counted skips are the triaged fixtures (dated,
//! decoded-reason debts) and a PER-ATOM |charge| > 1 (the
//! GetAtomChargeType guard; multi-charge totals — the multi-onium
//! /q+2..+4 and polyanion /q-2..-4 families — are compared since the
//! 2026-10-31 census). No silent skips; a triaged fixture that
//! passes is a retirement to act on.

use std::fs;
use std::path::PathBuf;

fn suite_dir() -> PathBuf {
    std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"))
}

#[test]
fn the_identifier_equals_the_oracle_recording() {
    let dir = suite_dir();
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
    let triage: serde_json::Value = match fs::read_to_string(dir.join("corpus/triage.json")) {
        Ok(text) => serde_json::from_str(&text).expect("valid triage"),
        Err(_) => serde_json::Value::Null,
    };

    let mut pass = 0;
    let mut skipped = 0;
    let mut retired = 0;
    let mut beyond = 0;
    let mut failures = Vec::new();
    let mut retirements = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        // Every fixture is compared; the triage file marks EXPECTED
        // failures (dated, decoded debts). A triaged fixture that
        // passes is a retirement to act on, never a silent skip.
        let expected_failure = triage.get(id).map(|_| ()).is_some();
        let recorded = case["inchi"].as_str().expect("inchi");

        let molfile = case["molfile"].as_str().expect("molfile");

        match inchi_rs::Engine::from_molfile(molfile).and_then(|engine| engine.inchi_string()) {
            Ok(assembled) => {
                // Full byte equality: every layer present, in the
                // oracle's order, byte-identical.
                let mut ok = assembled == recorded;
                let mut failures_now = Vec::new();
                if !ok {
                    failures_now.push(format!("assembled {assembled:?} != oracle {recorded:?}"));
                }
                // The native hasher must reproduce every recorded
                // InChIKey from the recorded identifier.
                let key = inchi_rs::key::inchikey(recorded).expect("keys");
                if key != case["inchikey"].as_str().unwrap_or("") {
                    ok = false;
                    failures_now.push(format!(
                        "key {key} != oracle {}",
                        case["inchikey"].as_str().unwrap_or("")
                    ));
                }
                if ok {
                    if expected_failure {
                        retired += 1;
                        retirements.push(id.to_string());
                    } else {
                        pass += 1;
                    }
                } else if !expected_failure {
                    for line in &failures_now {
                        failures.push(format!("{id}: {line}"));
                    }
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    beyond += 1;
                } else if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}\n(pass: {pass}, expected-failure skips: {skipped},          RETIRED (triaged, now passing — update the triage file): {retired}, \
         beyond-±1: {beyond}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "gate: pass {pass}, expected-failure skips {skipped}, retired {retired}, \
         beyond-±1 {beyond}, total {}",
        cases.len()
    );
    assert!(
        retirements.is_empty(),
        "{} triaged fixture(s) now PASS — retire their triage entries:\n{}",
        retirements.len(),
        retirements.join("\n")
    );
}
