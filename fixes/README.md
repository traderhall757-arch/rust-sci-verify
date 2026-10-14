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
| [puruspe/0001-Fix-betai-for-large-a-and-b.patch](puruspe/0001-Fix-betai-for-large-a-and-b.patch) | [001](../findings/001-puruspe-betai.md) | puruspe 0.4.4 = `main` at `f3eb94b` | revised (fix 3 now switches at the mean); awaiting the owner's approval to submit |

## Notes for submission (puruspe)

- Upstream CI runs `cargo fmt --check` and `cargo clippy -- -D warnings`. Both
  already fail on upstream `main` with Rust 1.97 (`src/beta.rs` formatting and
  `clippy::unnecessary_min_or_max` in `src/bessel.rs`), independently of this
  patch; the patch's own lines are formatted.
- The commit author is a placeholder (`Claude <noreply@anthropic.com>`); the
  owner would normally submit under their own name (`git am` then
  `git commit --amend --reset-author`).

## Draft pull request text (puruspe)

For the owner to post if approved; the commit message in the patch carries
the technical detail.

> **Fix `betai`: wrong values for large a and b and near the switch point**
>
> `betai` returns wrong values in three situations, each with a one-line
> cause: 0 or -0 whenever a and b both exceed 3000 (e.g.
> `betai(1e4, 1e4, 0.5)` is 0, should be 0.5), and values off by 2% up to
> a factor of 1e37 near the switch between the two continued fractions
> (e.g. `betai(1e3, 500.0, 0.6)` is 1.05, should be 4.8e-8;
> `betai(1.0, 0.5, 0.9999999999)` is 0.197, should be 0.99999). The commit
> message explains each cause and how the fix was checked against mpmath
> (20,000 random inputs plus fixed grids). A regression test is added to
> `tests/beta_test.rs`.
>
> Not changed: the general accuracy of about 1e-10 (up to about 1e-9 in
> some cases, and up to 5e-6 on results below 1e-115 with a, b > 3000),
> which comes from the `ln_gamma`-based prefactor.
>
> AI disclosure: the analysis and the patch were written by Claude (an AI
> model by Anthropic) for the rust-sci-verify project, which checks Rust
> numerical crates against mpmath; submitted by me on the project's behalf.

