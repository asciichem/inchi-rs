//! The gate: the same inchi-tests corpus the reference engine passes.
//! Contract: the engine's emitted segments equal the oracle's layer
//! prefix byte-for-byte — each milestone extends the prefix.
//! Mobile-H groups and the m3-tagged S/P slice gate live since M3
//! slice 2; counted skips are the triaged fixtures (dated,
//! decoded-reason debts) and isotopic fixtures (/i). |charge| > 1
//! counts beyond-M1. No silent skips.

use std::fs;
use std::path::PathBuf;

fn suite_dir() -> PathBuf {
    std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"))
}

#[test]
fn emitted_layers_agree_with_the_oracle_prefix() {
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
    let mut isotopic = 0;
    let mut beyond = 0;
    let mut failures = Vec::new();
    let mut retirements = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        // Every fixture is compared; the triage file marks EXPECTED
        // failures (dated, decoded debts). A triaged fixture that
        // passes is a retirement to act on, never a silent skip.
        let expected_failure = triage.get(id).map(|_| ()).is_some();
        let recorded: Vec<String> = case["inchi"]
            .as_str()
            .expect("inchi")
            .split('/')
            .map(str::to_string)
            .collect();
        let isotopic_scope = recorded.iter().any(|segment| segment.starts_with('i'));
        if isotopic_scope {
            // Isotopes ride their own milestone (/i layer); the reader
            // also sees explicit D atoms as heavy ones until then.
            isotopic += 1;
            continue;
        }

        let molfile = case["molfile"].as_str().expect("molfile");

        match inchi_rs::Engine::from_molfile(molfile).and_then(|engine| engine.segments()) {
            Ok(segments) => {
                // Kind-matched: each emitted segment is compared to
                // the oracle's recording of THAT layer kind, wherever
                // InChI orders it (formula = the bare first segment;
                // others by their letter).
                let mut ok = true;
                let mut failures_now = Vec::new();
                for (kind, segment) in &segments {
                    let expected = match kind {
                        None => recorded.get(1),
                        Some(letter) => recorded
                            .iter()
                            .find(|s| s.starts_with(*letter) && s.len() > 1),
                    };
                    match expected {
                        Some(expected) if expected == segment => {}
                        expected => {
                            ok = false;
                            failures_now.push(format!(
                                "{kind:?} segment {segment:?} != oracle {expected:?}"
                            ));
                        }
                    }
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
        "{} layer failure(s):\n{}\n(pass: {pass}, expected-failure skips: {skipped},          RETIRED (triaged, now passing — update the triage file): {retired}, isotopic skips: \
         {isotopic}, beyond-±1: {beyond}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "gate: pass {pass}, expected-failure skips {skipped}, retired {retired}, isotopic \
         {isotopic}, beyond-±1 {beyond}, total {}",
        cases.len()
    );
    assert!(
        retirements.is_empty(),
        "{} triaged fixture(s) now PASS — retire their triage entries:\n{}",
        retirements.len(),
        retirements.join("\n")
    );
}
