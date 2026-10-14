# Status

_Updated at the end of every session. Start here._

## Done

- **Milestone 1a** (PR from `claude/eloquent-einstein-7aly2r`): workspace, `harness`
  crate, `reference/` generator (mpmath 1.4.1 at 50+100 digits, SciPy 1.18.1,
  NumPy 2.5.4, Python 3.13, uv lockfile), data format schema 1, first dataset
  `erfcx` (361 points), `check-puruspe` (erfcx: 317 pass, 44 overflow
  convention, 0 wrong), CI jobs fmt / clippy / test / reference-fresh.
  Local: line coverage 99.2%, harness mutation score 100% (44/44 viable).

## Next

1. Owner merges the 1a PR.
2. **Milestone 1b** on the same branch, restarted from `main`: coverage job (90%
   lines), guardrail scripts (orphan files via dep-info, public API vs reference
   tests, cargo-machete, pin check), cargo-mutants on PR diffs with score in the
   job summary, weekly workflow (latest crates + SciPy in a separate environment,
   issue in this repo only), SessionStart hook for the Python environment.
   Then propose a minimum mutation score to the owner.
3. Milestone 2 in a fresh session.

## Open decisions (owner)

- After 1b: minimum mutation score for a required check.
- After 1b's first CI run: set up the `main` ruleset (steps in the 1a PR).

## Known problems / notes

- SciPy values may differ in the last bit across platforms; CI (ubuntu-24.04) is
  the source of truth for `reference/data/`.
- mpmath's `erfc` overflows internally for x > ~1e154; `erfcx` uses an enveloping
  asymptotic series for x >= 1e6 (tested against the direct formula).
