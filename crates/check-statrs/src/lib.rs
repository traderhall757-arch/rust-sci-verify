//! Verification of [statrs](https://crates.io/crates/statrs) against the
//! committed reference data. Everything lives in the test module; run with
//! `cargo test -p check-statrs`, which also rewrites
//! `reports/statrs-<version>.md`.

#[cfg(test)]
mod expected;
#[cfg(test)]
mod verify;
