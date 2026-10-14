# Status

_Updated at the end of every session. Start here._

## Done

- **Milestone 1a** ([PR #1](https://github.com/traderhall757-arch/rust-sci-verify/pull/1), merged):
  workspace, `harness` crate, `reference/` generator (mpmath 1.4.1 at 50+100
  digits, SciPy 1.18.1, NumPy 2.5.4, Python 3.13, uv lockfile), data format
  schema 1, dataset `erfcx` (361 points), `check-puruspe` (erfcx: 317 pass, 44
  overflow convention, 0 wrong), CI jobs fmt / clippy / test / reference-fresh.
- **Milestone 1b**: CI jobs `coverage` (99.2% lines, minimum 90%), `guardrails`
  (orphan files, API vs reference checks, cargo-machete, pins, script tests),
  `mutants` (PR diff: no surviving mutants unless excluded with a reason; weekly: at least 95%); weekly workflow (latest releases, issue in
  this repo only, full mutation run); SessionStart hook. See `docs/guardrails.md`.

## Next

1. Owner sets up the `main` ruleset requiring the CI checks.
2. **Milestone 2** in a fresh session: reproduce the four pilot findings
   (betai, gammp, statrs gamma pole, statrs erfc) plus the SciPy K_100 false
   alarm, with findings files. New datasets: betainc, gammainc, gamma, erfc, kv/yv.

## Open decisions (owner)

- None.

## Known problems / notes

- SciPy values may differ in the last bit across platforms; CI (ubuntu-24.04) is
  the source of truth for `reference/data/` (so far CI reproduced it exactly).
- mpmath's `erfc` overflows internally for x > ~1e154; `erfcx` uses an enveloping
  asymptotic series for x >= 1e6 (tested against the direct formula).
- The weekly job has been dry-run locally (no-difference and difference cases);
  its first real run will be the first Monday after merge, or run it by hand.
