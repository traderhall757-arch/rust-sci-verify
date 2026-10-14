# Guardrails

Every guardrail runs in CI (`.github/workflows/ci.yml`); a violation fails the
pull request. Each one is listed with how it works and why.

| # | Guardrail | CI job | How |
|---|---|---|---|
| 1 | No orphan `.rs` files | `guardrails` | `scripts/check_orphans.py` builds everything (`cargo build` and `cargo test --no-run`, all targets) and reads the compiler's own dep-info files for the artifacts of *that* run. Any `.rs` under `crates/` not listed is an orphan. |
| 2 | Warnings are errors | `clippy`, all jobs | `cargo clippy --all-targets -- -D warnings`, `RUSTFLAGS=-D warnings`, and workspace lints deny `dead_code`, `unused`, `unreachable_pub`, `missing_docs`, `unsafe_code`. |
| 3 | No stubs | `clippy` | Workspace lints deny `clippy::todo`, `clippy::unimplemented`, `clippy::dbg_macro`. |
| 4 | Coverage at least 90% of lines | `coverage` | `cargo llvm-cov --workspace --fail-under-lines 90`; the table goes to the job summary. |
| 5 | Mutation testing | `mutants` (PRs), weekly `mutants-full` | `cargo-mutants` on the PR's changed Rust code; all code weekly. Score in the job summary (`scripts/mutation_score.py`). Report-only until the owner sets a minimum. |
| 6 | Every public function has a reference test | `guardrails` | `scripts/check_api_tests.py`; see below. |
| 7 | No unused dependencies | `guardrails` | `cargo machete`. |
| 8 | Pinned versions | `guardrails` | `scripts/check_pins.py`: exact toolchain in `rust-toolchain.toml`; `Cargo.lock` and `reference/uv.lock` committed; Python deps `==`; target crates `=x.y.z`; CI tools `@version`. |

The guardrail scripts have their own tests (`scripts/tests/`), run in the
`guardrails` job.

## Coverage: what is counted, and exclusions

`cargo llvm-cov` counts every line compiled in the workspace, test modules
included. **There are no exclusions.** If one is ever needed (for example a
line that is reachable only when the disk is full), it must use
`#[coverage(off)]` or be listed here with the reason, and it appears in the PR's
"Checks changed" section.

## "Every public function has a reference test"

Each crate declares a role in `Cargo.toml` under `[package.metadata.rsv]`:

- **implementation**: scientific code we write (Phases B and C). Every `pub fn`
  must be named by a check (`function: "crate::path::name"`) in a check crate,
  i.e. verified against mpmath/SciPy reference data.
- **check**: verifies one target crate; it defines no `pub fn`.
- **harness**: test tooling; it contains no scientific functions. Its own
  correctness is guarded by coverage and mutation testing instead.

There are no implementation crates yet, so today this guardrail enforces the
roles and the "check crates export nothing" rule; it starts biting in Phase C.

## Weekly job (`.github/workflows/weekly.yml`)

Mondays 06:23 UTC, or by hand from the Actions tab ("Run workflow").

1. Lets every check crate resolve the latest release of its target crate
   (`scripts/relax_pins.py`, then `cargo update`), and upgrades SciPy, NumPy and
   mpmath in a separate Python environment.
2. Regenerates the reference data and runs `cargo test` in the runner's
   working copy. **Nothing is committed**: the job can only read the repository.
3. `scripts/weekly_compare.py` compares with the committed data and reports. If
   anything differs, it opens an issue labelled `weekly-verification` in this
   repository, or updates the open one. It never posts anywhere else.
4. Runs mutation testing on all code and reports the score.
