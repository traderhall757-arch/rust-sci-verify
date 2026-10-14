# Status

_Updated at the end of every session. Start here._

## Done

- **Milestone 1a** ([PR #1](https://github.com/traderhall757-arch/rust-sci-verify/pull/1), merged):
  workspace, `harness` crate, `reference/` generator (mpmath 1.4.1 at 50+100
  digits, SciPy 1.18.1, NumPy 2.5.4, Python 3.13, uv lockfile), data format
  schema 1, dataset `erfcx`, `check-puruspe`, CI jobs fmt / clippy / test /
  reference-fresh.
- **Milestone 1b** ([PR #2](https://github.com/traderhall757-arch/rust-sci-verify/pull/2), merged):
  CI jobs `coverage` (minimum 90% of lines), `guardrails`, `mutants` (PR diff: no
  surviving mutants unless excluded with a reason; weekly: at least 95%); weekly
  workflow; SessionStart hook. See `docs/guardrails.md`.
- **Milestone 2** (reproduce the pilot findings), one PR per function family:
  1. **erfc** ([PR #3](https://github.com/traderhall757-arch/rust-sci-verify/pull/3)): statrs `erfc`, finding 004 (low: about 10 digits; 111 of 281 points); new `check-statrs` crate; harness `Expected` lists all points of a finding.
  2. **gamma** ([PR #4](https://github.com/traderhall757-arch/rust-sci-verify/pull/4)): statrs `gamma`, finding 003 (high: finite at poles, early overflow; 187 of 901 points).
  3. **gammainc** ([PR #5](https://github.com/traderhall757-arch/rust-sci-verify/pull/5)): puruspe `gammp`, finding 002 (medium: NaN for a >= 100; 78 of 192 points); tolerance 1e-12. Also the CI fix: the `mutants` job passes when changed code has nothing to mutate.
  4. **betainc** ([PR #6](https://github.com/traderhall757-arch/rust-sci-verify/pull/6)): puruspe `betai`, finding 001 (high: 0 instead of 0.5 at a = b = 1e4; 336 of 639 points); tolerance 1e-12.
  5. **kv** (last PR): puruspe `Kn` is right where SciPy is wrong: K_100(0.0631) = inf in SciPy, 5.875e305 in mpmath and puruspe (the pilot quoted 5.92e305 at x = 0.063096); 8 "SciPy wrong" points, no finding.

  All four findings are **approved by owner** (2026-10-14); none is reported
  upstream (that needs a separate approval per report).

## Next

1. Owner merges the last Milestone 2 PR (kv). Milestone 2 is then complete; start a fresh session.
2. Owner sets up the `main` ruleset requiring the CI checks (if not done yet).
3. **Milestone 3** in a fresh session: diagnose `puruspe::betai` (root cause and a
   minimal, tested, fork-ready fix; not submitted). Leads are in finding 001.
4. Milestone 4 adds `yv` (puruspe Bessel functions over wider ranges), with its
   own dataset and checks (owner: fine as long as it is tested).

## Decisions recorded (owner, 2026-10-14)

- Findings 001 to 004 approved.
- Tolerance 1e-12 accepted for `gammainc` and `betainc` (1e-13 elsewhere);
  **revisit later** (owner's note), e.g. once a Boost/Temme-quality
  implementation shows what is achievable.
- Milestone 2 split into one PR per function family.

## Known problems / notes

- SciPy values may differ in the last bit across platforms; CI (ubuntu-24.04) is
  the source of truth for `reference/data/` (so far CI reproduced it exactly).
- mpmath's `erfc` overflows internally for x > ~1e154; `erfcx` and `erfc` use an
  enveloping asymptotic series for x >= 1e6 (tested against the direct formula).
- `mpmath.betainc` takes minutes per point near the mode for a, b around 1e4, so
  the `betainc` dataset uses its own positive-term series (DLMF 8.17.8), tested
  against `mpmath.betainc` and the exact binomial sum.
- Under the error metric, a subnormal result flushed to zero is wrong (relative
  to the smallest normal double). SciPy does this for erfc(26.75), erfc(27.0)
  and gammainc(200, 2).
- The weekly job has been dry-run locally; its first real run is the first
  Monday after merge, or run it by hand from the Actions tab.
