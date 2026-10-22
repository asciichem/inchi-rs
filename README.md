# inchi-rs

The production InChI engine: a native Rust implementation of InChI
derivation (graph canonicalization and layer serialization), with no
third-party binary at runtime.

It is gated by the shared oracle corpus
[inchi-tests](https://github.com/asciichem/inchi-tests), pinned in CI
at `v0.1.0`. The corpus is recorded from IUPAC `inchi-1` (MIT) at
development time; the IUPAC binary is a test oracle only and is never
shipped or executed by this engine.

Status: M0 (repo, reader, gate), M1 (formula + `/p`) and M2 (canonical
ranking + `/c`) green: the native canonicalization reproduces the
oracle's connection segments across the corpus, with the single
symmetric-system disagreement triaged in the corpus. M3 (mobile-H) is
next. See [ARCHITECTURE.md](ARCHITECTURE.md) for the module map and
extension axes.

## Commands

```sh
INCHI_TESTS=/path/to/inchi-tests cargo test   # unit tests + oracle gate
INCHI_TESTS=/path/to/inchi-tests cargo run --example refinement-census
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

MIT.
