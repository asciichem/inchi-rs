# inchi-rs

The production InChI engine: a native Rust implementation of InChI
derivation (graph canonicalization and layer serialization), with no
third-party binary at runtime.

It is gated by the shared oracle corpus
[inchi-tests](https://github.com/asciichem/inchi-tests), pinned in CI
at `v0.1.0`. The corpus is recorded from IUPAC `inchi-1` (MIT) at
development time; the IUPAC binary is a test oracle only and is never
shipped or executed by this engine.

Status: M0 (repo, reader, gate) and M1 (formula + `/p`) green over
both the identity seeds and the 1500-fixture breadth corpus; M2
(ranking and `/c`) in progress. See [ARCHITECTURE.md](ARCHITECTURE.md)
for the module map and extension axes.

## Commands

```sh
INCHI_TESTS=/path/to/inchi-tests cargo test   # unit tests + oracle gate
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

MIT.
