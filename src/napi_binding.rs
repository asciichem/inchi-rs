//! The Node binding (feature `napi`): the engine for Node consumers
//! via napi-rs.

use napi_derive::napi;

/// Derives the identifier from a V2000 molfile.
#[napi]
pub fn from_molfile(molfile: String) -> napi::Result<String> {
    crate::Engine::from_molfile(&molfile)
        .and_then(|engine| engine.inchi_string())
        .map_err(|error| napi::Error::from_reason(error.message().to_string()))
}

/// Derives the InChIKey from a V2000 molfile.
#[napi]
pub fn inchikey(molfile: String) -> napi::Result<String> {
    crate::Engine::from_molfile(&molfile)
        .and_then(|engine| engine.inchikey())
        .map_err(|error| napi::Error::from_reason(error.message().to_string()))
}
