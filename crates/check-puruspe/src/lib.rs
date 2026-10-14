//! Verification of [puruspe](https://crates.io/crates/puruspe) against the
//! committed reference data. Everything lives in the test module; run with
//! `cargo test -p check-puruspe`, which also rewrites
//! `reports/puruspe-<version>.md`.

#[cfg(test)]
mod verify;
