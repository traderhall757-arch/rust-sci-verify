//! Shared verification harness for rust-sci-verify.
//!
//! A *check crate* (one per target crate, e.g. `check-puruspe`) calls
//! [`Check::run`] once per function with a closure that evaluates the target
//! crate at a reference [`Point`], then passes all results to [`finish`].
//! `finish` writes the markdown report to `reports/` and panics if the set of
//! points where the crate is wrong differs from the declared
//! [`Expected`] failures, each of which must cite a file in `findings/`.
//!
//! Reference data lives in `reference/data/*.jsonl` (format:
//! `docs/data-format.md`); the classification rules are in
//! `docs/conventions.md`.

mod check;
mod compare;
mod data;
mod lockfile;
mod report;

use std::path::{Path, PathBuf};

pub use check::{Check, Expected, FunctionResult, Outcome, PointResult, expectation_problems};
pub use compare::{Verdict, classify, rel_error};
pub use data::{DataError, Dataset, Header, Point, Tag, Tolerance, ToolVersion, Values};
pub use lockfile::locked_version;
pub use report::render;

/// The workspace root (two levels above this crate's manifest).
#[must_use]
pub fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// Writes the report for `crate_name` and panics unless every crate-wrong
/// point is exactly matched by an [`Expected`] entry that cites an existing
/// findings file.
///
/// # Panics
///
/// On any unexpected failure, unexpected pass, missing findings file, or if
/// the report or `Cargo.lock` cannot be read or written.
pub fn finish(crate_name: &str, results: &[FunctionResult]) {
    let root = workspace_root();
    let lock = std::fs::read_to_string(root.join("Cargo.lock")).expect("read Cargo.lock");
    let version = locked_version(&lock, crate_name).expect("crate version in Cargo.lock");
    let path = root
        .join("reports")
        .join(format!("{crate_name}-{version}.md"));
    std::fs::create_dir_all(root.join("reports")).expect("create reports/");
    std::fs::write(&path, render(crate_name, &version, results)).expect("write report");
    let problems = expectation_problems(results, &root.join("findings"));
    assert!(
        problems.is_empty(),
        "{} problem(s) (full report: {}):\n{}",
        problems.len(),
        path.display(),
        problems.join("\n")
    );
}
