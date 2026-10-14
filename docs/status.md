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
- **Milestone 2 complete** ([PR #7](https://github.com/traderhall757-arch/rust-sci-verify/pull/7), merged).
- **Milestone 3** (diagnose `puruspe::betai`; [PR #8](https://github.com/traderhall757-arch/rust-sci-verify/pull/8), merged): root cause in
  [finding 001](../findings/001-puruspe-betai.md). Three bugs cause the wrong
  answers: missing parentheses in the quadrature integrand (every input with
  a, b > 3000 gives 0 or -0), a branch chosen by the sign of an underflowed
  value (-0 instead of 1), and a mistyped continued-fraction switch point
  (3% error at (1e4, 1e3, 0.9)). A fourth issue, the prefactor built from
  puruspe's 7-term `ln_gamma`, causes the 1e-11 to 6e-10 accuracy loss and is
  not fixed. New:
  - `fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch`: three-line fix plus a
    regression test, against puruspe 0.4.4 (= upstream `main`). **Not submitted.**
  - `crates/puruspe-fix`: verbatim copy of the patched code, checked as
    `puruspe_fix::betai`: no wrong answers left; worst error 5.15e-10 (CI bound
    6e-10). Report `reports/puruspe-fix-0.1.0.md`.
  - Dataset `betainc_large` (144 points, a, b in 3001..1e5): puruspe 0.4.4 is
    wrong at all of them.
  - Guardrail 9, `scripts/check_fixes.py`: the fix crate must equal the patch,
    and the patched crate's own tests must pass (CI `guardrails` job).
  - Mutation testing now runs the whole workspace's tests per mutant
    (`test_workspace = true`); 22 equivalent mutants in the copied upstream
    code are excluded with reasons in `.cargo/mutants.toml`.

- **Milestone 3 follow-up** (revise fix 3;
  [PR #10](https://github.com/traderhall757-arch/rust-sci-verify/pull/10), merged): a 20,000-point random check (SciPy, then mpmath for every
  difference above 1e-9) showed that puruspe 0.4.4's switch between its two
  continued fractions also gives **wrong answers at ordinary parameters**
  (a or b <= 3000): up to a factor of 1e37 when 1/a + 1/b < 1, up to 80% for
  b < 2 with x near 1 (e.g. `betai(1, 0.5, 0.9999999999)` = 0.197, true
  0.99999). The first fix 3 lost accuracy for small b (up to 1.1e-8).
  - Fix 3 now switches at the mean a / (a + b) (as Boost.Math does for
    a, b >= 1): no wrong answers left; worst 1.1e-9 on the random check
    where a or b <= 3000 (first fix 1.1e-8). Up to 4.8e-6 on results below
    1e-115 with a, b > 3000: the known prefactor floor (bug 4).
  - New dataset `betainc_switch` (186 points: small b with x = 1 - y, and
    points between the old and Numerical Recipes switch points). puruspe
    0.4.4 fails 170 (66 wrong answers); the fix fails 144, worst 3.0e-9.
  - Patch regenerated; its regression test now also fails with the first
    fix. The fix crate has per-dataset error bounds (2.5e-10 `betainc`,
    tightened from 6e-10; 5e-10 `betainc_large`; 3.5e-9 `betainc_switch`).
  - Draft upstream PR text with the disclosure line "submitted by me on the
    project's behalf" in `fixes/README.md`.

- **Milestone 4, part 1: `yv`** (Bessel Y; [PR #12](https://github.com/traderhall757-arch/rust-sci-verify/pull/12), merged):
  - Datasets `yv` (120 points, x <= 20, tolerance 1e-13) and `yv_large`
    (80 points, x from 50 to 1e6, tolerance 5e-10), integer orders 0..500.
    Tolerances come from the condition number |x Y'/Y| (tested); SciPy's
    worst is 5.0e-14 and 2.6e-11.
  - `puruspe::Yn` passes both (worst 3.4e-14 and 5.5e-11).
  - `puruspe::Jnu_Ynu` fails 29 points: **finding 005**, approved by owner (2026-10-14). It panics
    for x above about nu + 9832 (Numerical Recipes' continued fraction hits its
    10,000-iteration cap) and loses accuracy as x grows (1.4e-9 at
    Y_500(1e4); up to 4e-7 off the grid).

- **Milestone 4, part 2: `jv`** (Bessel J; [PR #14](https://github.com/traderhall757-arch/rust-sci-verify/pull/14), merged):
  - Datasets `jv` and `jv_large`: the `yv` grids and tolerances (1e-13 and
    5e-10, from the condition number |x J'/J|, tested). SciPy returns 0 for
    J_200(5) = 4.8e-296 ("SciPy wrong"); its worst on `jv_large` is 2.0e-10.
  - `puruspe::Jn` fails 27 points: **finding 006** (high), approved by
    Claude after the checklist. Wherever the true J_n(x) underflows (x well
    below n) it returns inf or an arbitrary number (J_500(1) = 0.0046,
    J_500(20) = 1.8e90). Root cause: puruspe's `ldexp` returns infinity
    instead of 0 on underflow; confirmed by swapping in a correct `ldexp`.
    Accurate where J is a normal double.
  - `puruspe::Jnu_Ynu`'s J result fails 55 points: finding 005 extended
    (approved by Claude after the checklist): the same panics, NaN where J
    underflows (19 points), accuracy loss up to 8.9e-13 on `jv`, 1.9e-8 at
    J_200(1e4), 1.3e-9 off the grid for non-integer orders at small x.
  - Correction to finding 005: the panic boundary is not exactly
    x = nu + 9845; panics start near nu + 9832 and alternate with values up
    to about nu + 9851.

- **Milestone 4, part 3: `iv`** (modified Bessel I; [PR #15](https://github.com/traderhall757-arch/rust-sci-verify/pull/15), open):
  - Datasets `iv` (the `yv` grid, 1e-13), `iv_large` (integer orders,
    x = 50..1000 across the overflow edge, 17 `overflow` tags, 2e-13) and
    `iv_nonint` (orders 0.25..200.75, x = 1e-5..700, 2e-13); tolerances from
    |x I'/I| (tested). SciPy returns inf for I_100(720) = 7.1e307.
  - `puruspe::In` fails 48 points: **finding 007** (high), approved by
    Claude after the checklist. Three bugs: inf or arbitrary numbers where
    I_n(x) is tiny (finding 006's `ldexp`, plus `ans` underflowing before
    normalisation), inf for every x > 709.78 (`exp(x)` overflows in `i0`),
    errors up to 1.2e-3 for x >> n (Miller start order too low).
  - `puruspe::Inu_Knu`'s I result fails 49 points: **finding 008** (high),
    approved by Claude after the checklist. NaN where I underflows; inf for
    x above about 714 (subnormal K in the Wronskian); accuracy loss from the
    accumulated recurrence factor `fact` (confirmed; finding 005's J loss
    has the same code) and from cancellation at small x for negative
    fractional order.

- **Milestone 4, part 4: standard normal distribution** (statrs `Normal`;
  [PR #16](https://github.com/traderhall757-arch/rust-sci-verify/pull/16), open):
  - Datasets `norm_cdf` and `norm_sf` (441 points each: x = -40..40 step 0.25
    plus decades, tolerance 4e-13 from the condition number |x phi/Phi|,
    which reaches 1408 at x = -37.5; tested) and `norm_ppf` (138 points, both
    tails, the bulk, 1/2 +- 10^-k; p = 0 and 1 tagged `pole`; 1e-13). The
    quantile reference is Newton's method on log Phi, tested against
    `sqrt(2) erfinv(2p - 1)` and Phi(ppf(p)) = p. SciPy flushes subnormal
    cdf values at x = -37.75 and -38 to 0 ("SciPy wrong").
  - `Normal::inverse_cdf` passes every point (worst 2.6e-16).
  - `Normal::cdf` and `Normal::sf` fail 155 points each (up to 1.03e-10):
    **finding 009** (low), approved by Claude after the checklist. Root cause
    is statrs's `erfc` (finding 004), now **confirmed**: rounding its 13
    interval constants to float, as Boost does, makes `erfc`, `cdf` and `sf`
    pass every grid point and a 20,000-point random scan.
  - Only mean 0, sd 1 is checked; other parameters use the same `erfc` call.

## Next

1. **Milestone 4, continued** (next session): more statrs distributions
   (CDF/SF/PPF vs `scipy.stats`), one family per PR. Suggested order:
   Student's t, chi-squared and gamma (they exercise statrs's incomplete
   gamma and beta, so expect overlap with findings 001 to 003 style issues),
   then the `special` crate, then statrs hypothesis tests. Optional small
   follow-up: check the K result of `Inu_Knu` on the `kv` grid (not yet
   checked; same code path as finding 008).
2. Nothing is submitted upstream until the whole project is done and tested
   (owner, 2026-10-14). The puruspe fork (`traderhall757-arch/puruspe`,
   branch `fix-betai`) still holds the old commit; replace it with the
   revised patch only when submission is approved.

## Decisions recorded (owner, 2026-10-14)

- Findings 001 to 004 approved; finding 005 approved (same day, after PR #12).
  Finding 006 and the J part of 005 approved by Claude (checklist); findings
  007, 008 and 009 likewise.
- **Claude approves future findings itself**, after the confirmation
  checklist in `findings/README.md` (reproduction, arbiter check, second
  opinion, off-grid check, root cause). Upstream contact still needs the
  owner's approval per item.
- Tolerance 1e-12 accepted for `gammainc` and `betainc` (1e-13 elsewhere);
  **revisit later** (owner's note), e.g. once a Boost/Temme-quality
  implementation shows what is achievable.
- Milestone 2 split into one PR per function family.
- Finding 001's extended scope (all a, b > 3000) confirmed.
- Submitting `fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch` to puruspe
  approved (this one patch only; any later report or patch needs its own approval).
  On hold after the spot-check; the revised patch (Milestone 3 follow-up)
  needs approval again.
- Finding 001's wider scope (wrong answers also for a or b <= 3000, near the
  switch point) confirmed.
- **No upstream submissions** (patches, issues, pull requests) until the
  project is finished and everything is tested; the revised puruspe patch
  waits until then and needs approval at that time.
- The owner merges PRs without reviewing the numerical details (unfamiliar
  with the field), so CI checks and the evidence in each PR carry the
  verification; keep them complete and honest.
- The `betai` accuracy floor (Boost-style prefactor, better `ln_gamma`; also
  helps `gammp`) becomes a later milestone, after Milestone 4.

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
- `yv` and `yv_large` take about 10 s together; `jv` and `jv_large` about 1 s; the three `iv` datasets about 1 s; the three `norm` datasets about 2 s. The `yv` grids leave out
  x = 0 (a pole; `Jnu_Ynu` documents a panic there).
- `betainc_large` takes about a minute to generate (mpmath series at a, b = 1e5);
  `betainc_switch` about 30 s. Its grid leaves out a = 1e5 with x = 1 - 1e-4,
  where the reference series needs millions of terms.
- Random spot-checks of `betai` are done with scratch scripts (not committed,
  rule 6); the method and seed are recorded in finding 001.
- Upstream puruspe `main` fails its own `cargo fmt --check` and `cargo clippy`
  on Rust 1.97, independently of our patch (`fixes/README.md`).
- The weekly job has been dry-run locally; its first real run is the first
  Monday after merge, or run it by hand from the Actions tab.
