//! Probe: run the engine on a molfile path given as arg.
use std::fs;

use chemomer::Engine;

fn main() {
    let path = std::env::args().nth(1).expect("molfile path");
    let input = fs::read_to_string(path).expect("read");
    let engine = Engine::from_molfile(&input).expect("parses");
    for (kind, segment) in engine.segments().expect("segments") {
        match kind {
            Some(letter) => println!("/{letter}{segment}"),
            None => println!("/{segment}"),
        }
    }
}
