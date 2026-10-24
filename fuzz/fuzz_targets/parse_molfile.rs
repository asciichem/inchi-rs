#![no_main]

//! Fuzz target: the V2000 reader must accept or reject arbitrary
//! input cleanly — no panics, no hangs. Run manually:
//!   cargo +nightly fuzz run parse_molfile

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = inchi_rs::parse_molfile(text);
    }
});
