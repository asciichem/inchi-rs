# inchi-rs — architecture

The Rust InChI engine. Design principles: **encapsulation** (every
concern behind a module boundary; callers see types, not column
arithmetic or rule tables), **open/closed** (new InChI layers,
invariant schemes, and input formats join by adding modules and
registering implementations — never by editing the engine core), and
**tests extracted from old InChI**: the oracle corpus in
[inchi-tests](https://github.com/asciichem/inchi-tests) is recorded
from the IUPAC reference implementation and gates this engine
byte-for-byte. The oracle is a development-time test oracle only —
nothing from the InChI Trust ships here.

## The derivation pipeline

```
molfile text ──▶ molfile::parse ──▶ graph::Structure ──▶ layers (registry) ──▶ segments
     (input edge)     (encapsulated)     (graph semantics)   (OCP axis)       (InChI order)
```

Each layer owns consecutive *decisions* of the derivation; the engine
only runs the registry and tags emitted segments by layer kind.
InChI's own segment order (formula, /c, /h, /q, /p, /i, /t…) is a
property of the format, not of the registry order — gates match
segments **by layer kind**, so a milestone-partial engine claims
exactly the layers it can decide.

## Module map

| Module | Responsibility | Extension axis |
|---|---|---|
| `error` | `Error` / `Result`; gate-visible messages. | — |
| `molfile` | V2000 subset reader → `Structure`. Column arithmetic and unsupported-construct reporting live here and nowhere else. | New input formats (SDF, SMILES) = new reader modules behind the same `Structure`. |
| `graph` | `Structure` / `AtomRecord`, bond orders, valences, implicit-H and composition semantics. Graph queries (bond-order sums, per-atom H) are graph concerns, not layer concerns. | New element/valence data lands in one table. |
| `layers` | One module per InChI layer; the `LayerEmitter` trait is the OCP seam: `name()`, `segment_kind()` (None = formula, Some('p')…), `emit(&Structure) -> Option<String>`. | A new layer = a new module + one registration in `default_layers()`. |
| `layers::standardization` | The ±1 proton rules decoded from the oracle (anion/cation proton accounting) — shared by the layers the rules couple (formula adjusts H, proton emits /p). One home, no duplication. | M3 replaces its arithmetic with full mobile-H rules. |
| `canonical` (M2) | Native port of the reference's base canonicalization path (MIT source read at development time, never compiled or called): invariants → partition → equitable refinement → two-pass individualization search → connection table → `/c` DFS emission. | Invariants and the table comparison are the corpus-arbitrated seams; M3 adds t-group nodes to the same partition machinery. |
| `engine` | The façade: `Engine::from_molfile`, `segments()` — runs the registry, tags segments by kind. Closed for modification. | — |

## Encapsulation boundaries

- **Callers never see molfile column arithmetic** — `molfile::parse`
  is the only code that indexes a molfile line.
- **Layers never see each other.** The formula layer cannot consult
  the proton layer; coupled rules live in `standardization` and are
  consumed as plain values.
- **The engine never branches on layer identity.** It runs the
  registry and tags results; `segment_kind` carries the identity.
- Rule tables (valences, proton standardization) are private to the
  module that owns the semantics.

## Tests extracted from old InChI

The suite (inchi-tests) is the extracted behaviour of the old
InChI — every fixture is `(molfile, InChI, InChIKey)` recorded from
the IUPAC reference binary:

1. **Gate (`tests/gate.rs`)** — the whole-suite contract: for every
   fixture, every emitted segment must match the oracle's recording
   of that layer kind, byte-for-byte. M3-scoped (m3-tagged) and
   triaged fixtures are counted skips with reasons; |charge| > 1
   counts beyond-M1. No silent skips.
2. **Unit tests** — per module: molfile offsets, valence/H inference,
   Hill ordering, proton-rule cases (ammonium vs quaternary N vs
   carbanion), each decoded rule pinned by the fixture that
   discovered it.
3. **Extraction at scale** — `inchi-tests/scripts/import-pubchem-breadth.rb`
   pulls PubChem structures through the gem's molfile writer, verifies
   the oracle, and grows the corpus; disagreements become triage
   entries or bug fixes (this loop found four real identity bugs in
   its first run).

## Milestones → modules

| Milestone | Modules | Gate |
|---|---|---|
| M0 corpus | inchi-tests | coverage review |
| M1 formula + proton | `layers::formula`, `layers::proton`, `layers::standardization` | kind-matched byte agreement, 100% of in-scope |
| M2 /c connection | `canonical` (invariants, refinement, ranking), `layers::connection` | ≥99.9% full-string on slice; disagreements triaged |
| M3 /h /q /p mobile-H | `layers::mobile_h`, standardization v2 | ≥99.9% on slice |
| M4 stereo /t /b /m /s | `layers::stereo` | ≥99.9% incl. stereo slices |
| M5 InChIKey | `key` (SHA-256 construction) | 100% wherever the string agrees |
| M6 bindings | PyO3 / napi / wasm crates wrapping this one | each runtime passes the same suite |
| M7 engine claims | engine registration + docs | full-corpus refresh green |
