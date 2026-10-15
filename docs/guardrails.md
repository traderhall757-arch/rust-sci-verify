# Guardrails

Every guardrail runs in CI (`.github/workflows/ci.yml`); a violation fails the
pull request. Each one is listed with how it works and why.

| # | Guardrail | CI job | How |
|---|---|---|---|
| 1 | No orphan `.rs` files | `guardrails` | `scripts/check_orphans.py` builds everything (`cargo build` and `cargo test --no-run`, all targets) and reads the compiler's own dep-info files for the artifacts of *that* run. Any `.rs` under `crates/` not listed is an orphan. |
| 2 | Warnings are errors | `clippy`, all jobs | `cargo clippy --all-targets -- -D warnings`, `RUSTFLAGS=-D warnings`, and workspace lints deny `dead_code`, `unused`, `unreachable_pub`, `missing_docs`, `unsafe_code`. |
| 3 | No stubs | `clippy` | Workspace lints deny `clippy::todo`, `clippy::unimplemented`, `clippy::dbg_macro`. |
| 4 | Coverage at least 90% of lines | `coverage` | `cargo llvm-cov --workspace --fail-under-lines 90`; the table goes to the job summary. |
| 5 | Mutation testing | `mutants` (PRs), weekly `mutants-full` | `cargo-mutants` on the PR's changed Rust code; all code weekly. Score in the job summary (`scripts/mutation_score.py`). Rule below. |
| 6 | Every public function has a reference test | `guardrails` | `scripts/check_api_tests.py`; see below. |
| 7 | No unused dependencies | `guardrails` | `cargo machete`. |
| 8 | Pinned versions | `guardrails` | `scripts/check_pins.py`: exact toolchain in `rust-toolchain.toml`; `Cargo.lock` and `reference/uv.lock` committed; Python deps `==`; target crates `=x.y.z`; CI tools `@version`. |
| 9 | A tested fix is the patch offered upstream | `guardrails` | `scripts/check_fixes.py --run-tests`: applies each `fixes/` patch to the pinned upstream source, checks that the items listed under `[package.metadata.rsv.fix]` in the fix crate equal the patched ones character for character, and runs the patched crate's own tests. See `fixes/README.md`. |
| 10 | Committed data equals a fresh regeneration | `reference-fresh`, weekly `reference-fresh` | `python -m rsv_reference` regenerates the data, then `git status` must show no change in `reference/data/` (untracked files included). Pushes to `main` and the weekly run regenerate every dataset; pull requests only the datasets `scripts/select_datasets.py` says they can have changed. See below. |

The guardrail scripts have their own tests (`scripts/tests/`), run in the
`guardrails` job.

## Mutation-testing rule (approved by the owner)

A *mutant* is a small deliberate bug inserted into our code; a good test suite
fails on it ("caught"). Mutants that do not compile ("unviable") are ignored.

- **Pull requests:** the `mutants` job fails if any mutant in the changed code
  survives (missed or timed out). The only way around it is to exclude that
  mutant in `.cargo/mutants.toml` with a `# reason:` comment, typically because
  the mutant does not change behaviour so no test could catch it. The
  `guardrails` job checks that every exclusion has a reason, and every change
  to the file is a "Checks changed" item.
- **Weekly full run:** fails if the score over all code is below 95%.
- **Which tests run:** every mutant runs the whole workspace's tests
  (`test_workspace = true`), because implementation crates such as
  `puruspe-fix` are tested by the reference checks in the check crates.

A fixed percentage is not used for pull requests because they are small: one
survivor out of five mutants would already be 80%.

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

## Which datasets a pull request regenerates

Regenerating everything takes over ten minutes and grows with every dataset,
so on a pull request `scripts/select_datasets.py` decides, and prints one line
per dataset saying whether it is regenerated and why (also in the job summary).
It is safe because `python -m rsv_reference` generates each dataset in its own
fresh process, so a dataset's file depends only on:

1. the Python environment and the generator: any change under `reference/`
   other than `functions.py`, `data/` and `tests/`, or to `ci.yml` or the
   script itself, regenerates **everything**;
2. the parts of `functions.py` that run in every process at import (imports,
   module-level statements, decorated functions, assignments that call
   anything but `Spec`, and everything those reach): any change regenerates
   **everything**;
3. the definitions its `Spec` reaches through names, transitively: if any of
   them changed (compared as syntax trees, so comments and formatting do not
   count), or its data file changed, or it is new, **that dataset** is
   regenerated.

Anything the analysis cannot follow (`globals()`, `eval`, a name bound twice,
an import outside a short list, a module that `functions.py` imports importing
it back) regenerates everything. The script fails the job if the datasets it
reads from the source differ from the imported registry, or if a data file has
no dataset (or a dataset no file), so nothing is skipped silently. Every
commit on `main` finishes a full regeneration (CI runs on `main` are not
cancelled), and the weekly workflow repeats it with the pinned versions and
opens a `reference-stale` issue if anything differs.

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
