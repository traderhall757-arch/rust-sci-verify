# Roadmap

Mission: strengthen Rust's scientific-computing ecosystem through verification.
mpmath at high precision is the arbiter; SciPy is the practical reference.
Confirmed bugs become reproducible reports and, when small, fixes. Later the
same harness verifies implementations we write.

Background survey (October 2026, settled; do not redo): repository
`traderhall757-arch/claude`, branch `claude/vibrant-hopper-iq8rj5`,
`docs/survey-2026-10.md` and `pilot/`.

## Phase A: the verification harness (now)

Design (approved by the owner, with amendments):

1. Rust workspace: `crates/harness` + one `crates/check-<crate>` per target crate.
2. `reference/` Python package with exact pins (SciPy, NumPy, mpmath) writes
   `reference/data/*.jsonl` with both mpmath and SciPy values per point.
3. Rust tests read committed data; `cargo test` needs no Python.
4. CI regenerates the data and fails if it differs.
5. Every point classified: pass / crate wrong / SciPy wrong / both wrong /
   convention (see conventions.md). Convention only for accepted values (NaN,
   infinity, error) at rule-tagged points; finite values there are crate wrong.
6. Tolerance per function, with a documented reason.
7. Generated markdown report per crate version in `reports/`.
8. `findings/NNN-crate-function.md` per confirmed bug: minimal Rust repro,
   expected values (mpmath and SciPy), affected ranges, severity, status
   (draft / approved by owner / reported / fixed upstream).

Milestones:

1. Scaffold, split into two PRs:
   - **1a**: workspace, harness, data format, fmt/clippy/test CI, data freshness.
   - **1b**: coverage, guardrails, mutation testing, weekly workflow, SessionStart
     hook. After 1b, propose a minimum mutation score to make it a required check.
2. Reproduce the four pilot findings in the harness, with findings files:
   puruspe `betai` (large a or b), puruspe `gammp(500, 100)` = NaN,
   statrs `gamma(-6.0)` = 5.94e12, statrs `erfc(0.5)` about 10 digits.
   Include the SciPy K_100(0.0631) = inf false alarm as a "SciPy wrong" point.
3. Diagnose `puruspe::betai`: root cause and a minimal, tested, fork-ready fix.
   Do not submit it.
4. Expand coverage: statrs distribution CDF/PPF/SF; the `special` crate; puruspe
   Bessel functions over wider ranges; then statrs hypothesis tests vs `scipy.stats`.
5. Gap map: commonly used SciPy functions without a trustworthy Rust
   implementation. Feeds Phases B and C.

Phase A keeps running on schedule after B and C start.

## Phase C: missing special functions (after A, unless the gap map says otherwise)

Fresnel, exponential integrals, sine/cosine integrals, Struve, zeta, later
hypergeometric functions. Scalar functions, so they use the Phase A harness
unchanged. The gap map picks the order. Natural upstream home: puruspe (issue
#20 has requested Fresnel integrals since 2024; the betai fix builds the
relationship). Port only from clearly licensed sources (Boost, published
papers); Cephes licensing is unclear.

## Phase B: spectral analysis

`welch`, `periodogram`, `spectrogram`, `csd`, `coherence` on rustfft. Needs
array comparison: the data format already stores every value as a list.
Before building, verify sci-rs's existing filters (butter, sosfilt, sosfiltfilt,
savgol). Natural upstream home: sci-rs (active, slow to release).

Anything built in B or C passes the same harness and guardrails.

## Automation

- CI on every push/PR: fmt, clippy, tests, data freshness, all guardrails.
- Weekly: latest target crates and latest SciPy in a separate environment;
  opens or updates an issue in this repository if results differ. Never
  modifies committed data; never posts to other repositories. Full mutation run.
- Cloud sessions: SessionStart hook installs the pinned Python environment.

## Guardrails (CI fails on violation)

1. Orphan `.rs` files (compared with compiler dep-info).
2. Warnings are errors; deny dead code, unused, unreachable pub.
3. No stubs: deny `clippy::todo`, `clippy::unimplemented`.
4. Coverage at least 90% lines (`cargo llvm-cov`), exclusions explained.
5. Mutation testing (`cargo-mutants`) on changed code in PRs; all code weekly.
6. Every public function has a reference test (script compares API vs tests).
7. No unused dependencies (`cargo-machete`).
8. Pinned versions: `rust-toolchain.toml`, Python lockfile, `Cargo.lock`.
