//! Refinement census over the oracle corpus: how often does the
//! equitable refinement already produce a total order (no
//! individualization needed), and how large are the symmetric cells
//! that remain?
//!
//! Run with INCHI_TESTS pointing at an inchi-tests checkout:
//!   cargo run --example refinement-census

use std::fs;
use std::path::PathBuf;

use chemomer::canonical::equitable_partition;

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
    let mut largest_cell = 0usize;
    let mut hist = std::collections::BTreeMap::<usize, usize>::new();

    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        for case in parsed {
            let molfile = case["molfile"].as_str().expect("molfile");
            let structure = chemomer::parse_molfile(molfile).expect("parses");
            let partition = equitable_partition(&structure);
            cases += 1;
            atoms += partition.at_number().len();
            *hist.entry(partition.cell_count()).or_insert(0) += 1;
            if partition.is_discrete() {
                discrete += 1;
            } else {
                symmetric += 1;
                let n = partition.at_number().len();
                let mut sizes = vec![0usize; n + 1];
                for position in 0..n {
                    sizes[partition.rank_of(partition.at_number()[position]) as usize] += 1;
                }
                largest_cell = largest_cell.max(sizes.into_iter().max().unwrap_or(0));
            }
        }
    }

    println!("fixtures: {cases} ({atoms} heavy atoms)");
    println!("discrete: {discrete}");
    println!("symmetric: {symmetric} (largest non-singleton cell: {largest_cell})");
    println!("cells-per-fixture histogram: {hist:?}");
}
