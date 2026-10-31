//! The identity pipeline benchmark: molfile parse → identifier →
//! key, over representative corpus shapes. Development-time
//! regression check (run: cargo bench), not a CI gate.

use criterion::{criterion_group, criterion_main, Criterion};

const WATER: &str = "water\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
const ALANINE: &str = "l-alanine\n  t\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  0  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  2  0  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n";

fn bench_identity(c: &mut Criterion) {
    let mut group = c.benchmark_group("identity");
    group.throughput(criterion::Throughput::Elements(1));
    for (name, molfile) in [("water", WATER), ("l-alanine", ALANINE)] {
        group.bench_function(format!("{name} inchi_string"), |b| {
            b.iter(|| {
                chemomer::Engine::from_molfile(molfile)
                    .and_then(|engine| engine.inchi_string())
                    .expect("derives")
            })
        });
        group.bench_function(format!("{name} inchikey"), |b| {
            b.iter(|| {
                chemomer::Engine::from_molfile(molfile)
                    .and_then(|engine| engine.inchikey())
                    .expect("derives")
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_identity);
criterion_main!(benches);
