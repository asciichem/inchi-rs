//! The Python binding (feature `pyo3`): the engine as a native
//! extension module `inchi_rs` with `from_molfile(molfile) ->
//! (inchi, inchikey)`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Derives `(inchi, inchikey)` from a V2000 molfile.
#[pyfunction]
fn from_molfile(molfile: &str) -> PyResult<(String, String)> {
    let engine = crate::Engine::from_molfile(molfile)
        .map_err(|error| PyValueError::new_err(error.message().to_string()))?;
    let inchi = engine
        .inchi_string()
        .map_err(|error| PyValueError::new_err(error.message().to_string()))?;
    let inchikey = engine
        .inchikey()
        .map_err(|error| PyValueError::new_err(error.message().to_string()))?;
    Ok((inchi, inchikey))
}

#[pymodule]
fn inchi_rs(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(from_molfile, module)?)?;
    Ok(())
}
