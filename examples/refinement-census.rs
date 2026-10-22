//! Refinement census over the oracle corpus: how often does invariant
//! refinement already produce a total order (no individualization
//! needed), and how large are the symmetric classes that remain?
//!
//! Run with INCHI_TESTS pointing at an inchi-tests checkout:
//!   cargo run --example refinement-census

use std::fs;
use std::path::PathBuf;

use inchi_rs::canonical::{refine, DefaultInvariants};

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let fixtures = dir.join("corpus/fixtures");
    let mut files: Vec<PathBuf> = fs::read_dir(&fixtures)
        .expect("corpus fixtures dir")
        .filter_map(|entry| {
            let path = entry.expect("dir entry").path();
            let name = path.file_name()?.to_string_lossy().to_string();
            (name.starts_with("identity") && name.ends_with(".json")).then_some(path)
        })
        .collect();
    files.sort();

    let mut cases = 0usize;
    let mut atoms = 0usize;
    let mut discrete = 0usize;
    let mut symmetric = 0usize;
    let mut largest_class = 0usize;
    let mut hist = std::collections::BTreeMap::<usize, usize>::new();

    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        for case in parsed {
            let molfile = case["molfile"].as_str().expect("molfile");
            let structure = inchi_rs::parse_molfile(molfile).expect("parses");
            let classes = refine(&structure, &DefaultInvariants);
            cases += 1;
            atoms += classes.ids().len();
            *hist.entry(classes.count()).or_insert(0) += 1;
            if classes.is_discrete() {
                discrete += 1;
            } else {
                symmetric += 1;
                let mut sizes = vec![0usize; classes.count()];
                for &id in classes.ids() {
                    sizes[id] += 1;
                }
                largest_class = largest_class.max(sizes.into_iter().max().unwrap_or(0));
            }
        }
    }

    println!("fixtures: {cases} ({atoms} heavy atoms)");
    println!("discrete: {discrete}");
    println!("symmetric: {symmetric} (largest non-singleton class: {largest_class})");
    println!("classes-per-fixture histogram: {hist:?}");
}
