# inchi-rs

The production InChI engine: a native Rust implementation of InChI
derivation (graph canonicalization and layer serialization), with no
third-party binary at runtime.

It is gated by the shared oracle corpus
[inchi-tests](https://github.com/asciichem/inchi-tests), pinned in CI
at `v0.1.0`. The corpus is recorded from IUPAC `inchi-1` (MIT) at
development time; the IUPAC binary is a test oracle only and is never
shipped or executed by this engine.

Status: M0-M2 green (formula, `/p`, canonical ranking + `/c`), and M3
slice 1 landed: the `/h` layer emits ±1-standardized attached
hydrogens over the tautomer-free corpus slice, with a corpus-derived
minimum-hydrogen tie-break at table-equal terminals. Remaining M3
scope: mobile-H groups (`(H,x,y)`) and the orbit-representative
disagreements triaged in the corpus (v0.1.2). See
[ARCHITECTURE.md](ARCHITECTURE.md) for the module map and extension
axes.

## Commands

```sh
INCHI_TESTS=/path/to/inchi-tests cargo test   # unit tests + oracle gate
INCHI_TESTS=/path/to/inchi-tests cargo run --example refinement-census
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

MIT.
