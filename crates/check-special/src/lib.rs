//! Verification of [special](https://crates.io/crates/special) against the
//! committed reference data. Everything lives in the test module; run with
//! `cargo test -p check-special`, which also rewrites
//! `reports/special-<version>.md`.

#[cfg(test)]
mod expected;
#[cfg(test)]
mod verify;
