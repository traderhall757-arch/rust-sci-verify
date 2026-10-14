# Fixes

Fixes for bugs in target crates, ready to offer upstream. **None is submitted**
unless the owner approves that specific submission (CLAUDE.md rule 8).

Each fix has two parts that CI keeps identical:

- `fixes/<crate>/NNNN-*.patch`: a `git format-patch` against the pinned
  release, with a regression test in the upstream crate's own style. Apply it
  with `git am` in a clone of the upstream repository.
- `crates/<crate>-fix`: the patched items copied verbatim, verified by the
  harness against the reference data like any target crate (report in
  `reports/<crate>-fix-<version>.md`).

`scripts/check_fixes.py --run-tests` (CI `guardrails` job) applies the patch to
the pinned upstream source, checks that the copied items match it character
for character, and runs the patched crate's own test suite.

| Fix | Finding | Upstream base | Status |
|---|---|---|---|
| [puruspe/0001-Fix-betai-for-large-a-and-b.patch](puruspe/0001-Fix-betai-for-large-a-and-b.patch) | [001](../findings/001-puruspe-betai.md) | puruspe 0.4.4 = `main` at `f3eb94b` | on hold: a random spot-check found fix 3 makes some results worse (finding 001); to be revised and re-approved |

## Notes for submission (puruspe)

- Upstream CI runs `cargo fmt --check` and `cargo clippy -- -D warnings`. Both
  already fail on upstream `main` with Rust 1.97 (`src/beta.rs` formatting and
  `clippy::unnecessary_min_or_max` in `src/bessel.rs`), independently of this
  patch; the patch's own lines are formatted.
- The commit author is a placeholder (`Claude <noreply@anthropic.com>`); the
  owner would normally submit under their own name (`git am` then
  `git commit --amend --reset-author`).
