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

- **Milestone 4, part 3: `iv`** (modified Bessel I; [PR #15](https://github.com/traderhall757-arch/rust-sci-verify/pull/15), merged):
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
  [PR #16](https://github.com/traderhall757-arch/rust-sci-verify/pull/16), merged):
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

- **Milestone 4, part 5: Student's t distribution** (statrs `StudentsT`;
  [PR #17](https://github.com/traderhall757-arch/rust-sci-verify/pull/17), merged):
  - Datasets `t_cdf` and `t_sf` (1007 points each: 19 df from 0.1 to 1e8 x
    53 x from -1e300 to 1e300) and `t_ppf` (513 points: the same df x 27 p
    from 5e-324 to 1 - 1e-15; p = 0, 1 tagged `pole`, 7 `overflow`). All
    tolerances 2e-13: cdf/sf from the condition number |x f/F| (worst 901 at
    df = 1e8, x = -30), ppf from |log|x|| (worst 689.6); both tested. The
    reference is a new continued fraction for I_x(a, 1/2) (`_mp_betainc_cf`,
    DLMF 8.17.22), because `mpmath.betainc` fails to converge for df >= 1e5;
    tested against the series, `mpmath.betainc`, closed forms for df = 1, 2,
    3, and the 2F1 form at 300 digits. SciPy is wrong at 5 cdf, 5 sf and 36
    ppf points (e.g. it returns 0 at p = 1/2 +- 1e-10 for df = 4).
  - `StudentsT::cdf` and `sf` fail 144 points each: **finding 010** (high),
    approved by Claude after the checklist. 0 for |x| > 1.34e154 (x^2
    overflows), exactly 0.5 for x^2/df below about 1e-15 (1 - h formed from a
    rounded h, and a 1e-15 cut-off in `beta_reg`), and errors growing with df
    (8.9e-8 at df = 1e8 on the grid, 4.3e-4 off it). All causes confirmed
    with a switchable copy of `beta_reg`.
  - `StudentsT::inverse_cdf` fails 219 points: **finding 011** (high),
    approved by Claude after the checklist. It **never returns** at 23 points
    (every p for df = 1e8: the Newton step's density overflows to inf and
    AS 109's step-shrinking loop has no limit), panics for p < 7.4e-163,
    cuts tails off near 1e8 (-9.2e7 for -3.2e19), returns sqrt(df 2^-53)
    near p = 1/2, and is off by up to 67% for large df (AS R83's absolute
    stopping rule).
  - Harness: new `Check::run_with_time_limit` and `Outcome::Timeout` ("no
    result within N s", always crate wrong; `docs/conventions.md`), so an
    endless loop in a crate cannot hang CI. `t_ppf` uses a 5 s limit.

- **Milestone 4, part 6: gamma and chi-squared distributions** (statrs
  `Gamma`, `ChiSquared`;
  [PR #18](https://github.com/traderhall757-arch/rust-sci-verify/pull/18), merged):
  - Datasets `gamma_cdf` and `gamma_sf` (538 points each: 16 shapes from 0.01
    to 1e6 x x from 0 to 1e300 plus a + k sqrt(a); tolerance 1e-12, the
    `gammainc` tolerance, since with rate 1 this is P(a, x)) and `gamma_ppf`
    (432 points: same shapes x 27 p; p = 1 tagged `pole`, p = 0 untagged
    because the quantile is 0; 2e-13 from the largest |log x| 691.3, tested).
    `mpmath.gammainc` does not converge at a = 1e6, so the reference is the
    series DLMF 8.7.1 / continued fraction DLMF 8.9.2, tested against
    `mpmath.gammainc`, the finite sum DLMF 8.4.10 up to a = 1e6, and closed
    forms. SciPy is wrong at 3 cdf/sf points (2.1e-10 at a = 1e6) and 6 ppf
    points (p = 5e-324 for a >= 100).
  - `ChiSquared(k)` is statrs's `Gamma(k/2, 1/2)`, so it is checked on the
    same data (df = 2a, at 2x, half its quantile): cdf/sf bit-identical to
    `Gamma`, quantile failing at partly different points.
  - `cdf`/`sf` fail 22 and 21 points: **finding 012** (low), approved by
    Claude after the checklist. Digits lost for shape >= about 750 (1.2e-9 on
    the grid, 3.4e-9 off it): the log prefactor a ln x - x - ln Gamma(a)
    cancels between terms near 1.3e7 (statrs's `ln_gamma` itself is within
    about 1 ulp); subnormal results are flushed to 0 by an `ax < -709.78`
    shortcut.
  - `inverse_cdf` fails 170 points for each distribution: **finding 013**
    (high), approved by Claude after the checklist. About 9.3e-10 for every
    quantile below 1e-9 (the convergence test is effectively absolute 1e-9),
    0.2% off at p = 1 - 1e-15 (solves cdf = p instead of sf = 1 - p), up to
    8% off in the far lower tail for large shape (100-step cap returns an
    unconverged iterate), plus finding 012's errors. No hangs; the checks
    still use `run_with_time_limit` and run before `t_ppf`.

- **Milestone 4, part 7: the `special` crate, incomplete beta** (special
  0.14.2 `Beta::inc_beta`, `Beta::inv_inc_beta`;
  [PR #19](https://github.com/traderhall757-arch/rust-sci-verify/pull/19), merged):
  - special is MIT OR Apache-2.0. It exports `Beta` (inc_beta, inv_inc_beta,
    ln_beta), `Gamma` (gamma, digamma, trigamma, inc_gamma, ln_gamma), and
    with default features `Elliptic` and `LambertW`. Its default build is
    `no_std`: `ln`, `exp`, `lgamma`, `tgamma` come from the libm crate. New
    crate `check-special` pins `special = "=0.14.2"` with default features.
  - New dataset `betaincinv` (1620 points: the betainc a, b x 20 p from 0 to
    1, including 5e-324 and the doubles nearest 1 - 1e-15; 2e-13 from the
    largest |log x| 700.0, tested). Reference: Newton on log I in
    s = log(x / (1 - x)) with the continued fraction, tested against the
    series and closed forms for b = 1, a = 1, a = b = 1/2. p = 0 and 1 give
    the support's ends (untagged). One point (a = 1, b = 2, p = 5e-324) is an
    exact rounding tie and is `unstable-reference`. SciPy fails 105 points
    (NaN, 2.2e-308 for underflowing quantiles, up to 911x off at
    p = 5e-324).
  - `inc_beta` on the three betainc datasets fails 182 of 969 points (worst
    4.0e-9): **finding 014** (low), approved by Claude after the checklist.
    Prefactor from separately rounded logs up to 2.2e6 (139 points), series
    cut-off too early when terms shrink slowly (12; up to 14.8 million terms,
    0.13 s per call), rounding over long series (lead, 31). Off the grid:
    2,362 of 10,000, all with max(a, b) >= 283.8, worst 6.0e-10.
  - `inv_inc_beta` fails 462 of 1620: **finding 015** (high), approved by
    Claude after the checklist. p < 1.5e-162: NaN start (p * p underflows),
    then a `debug_assert!` panic in debug builds (CI) and an **endless
    loop** in release builds (161 points). 7 more endless loops in the far
    tails of large parameters (Newton step -inf, no iteration limit).
    Quantiles below about 4.6e-14 wrong (109 points, up to a factor of
    1e188: absolute stopping rule, start clamped to 1e-4, Newton in x),
    16 more off by up to 1.9e-4, 169 small (mostly finding 014's error).
    Off the grid 1,520 of 3,000 fail, 809 never return. It is the original
    of statrs's `inv_beta_reg` (finding 011: same causes).

- **Milestone 4, part 8: the `special` crate, gamma family** (special
  0.14.2 `Gamma::gamma`, `ln_gamma`, `inc_gamma`;
  [PR #20](https://github.com/traderhall757-arch/rust-sci-verify/pull/20), merged):
  - New datasets `gammaln` (956 points: the `gamma` grid, decades to 1e300,
    the overflow edge near 2.55e305, 1 +- d and 2 +- d for d = 1e-5..1e-15,
    large negative non-integers; 2e-13 from the reflection formula's
    cancellation, 664x at x = -2.75; tested), `gammasgn` (its sign, the
    `gammaln` grid without poles, tolerance 0) and `gammainc_large` (P(a, x)
    for a = 1e8..1e20, 51 points, 1e-12; tested against numerical
    integration). SciPy is wrong next to 1 and 2 for `gammaln` (up to 29%)
    and in the lower tail of `gammainc_large` (up to 53%).
  - `gamma` (libm `tgamma`), `ln_gamma` and its sign (libm `lgamma_r`)
    pass every point.
  - `inc_gamma` fails 46 of 781 points: **finding 016** (high), approved by
    Claude after the checklist. It returns exactly 1 for every x > 1e8
    (AS 239's shortcut, kept after special removed the p > 1000 branch that
    preceded it), so it is wrong for every shape above about 1e8; digits
    lost for large shape (log prefix from rounded terms; from a = 800 off
    the grid, 1.2e-9 at 1e6, 2.0e-7 at 1e8); subnormal results flushed to 0.
    All three causes confirmed with a switchable copy.

- **Milestone 4, part 9: the `special` crate, digamma and trigamma** (special
  0.14.2 `Gamma::digamma`, `Gamma::trigamma`;
  [PR #21](https://github.com/traderhall757-arch/rust-sci-verify/pull/21), merged):
  - Harness: new `Check::run_isolated` and `Outcome::Crash` ("crash: stack
    overflow", always crate wrong; `docs/conventions.md`, "Calls that end the
    process"). A stack overflow aborts the process and cannot be caught, so
    such a check evaluates its points in child processes (the same test
    executable running only that check); a dead child's point is recorded and
    a new child continues at the next point. The check goes first in its test.
  - New datasets `digamma` (949 points: the `gammaln` grid's ranges plus the
    double nearest the zero x0 = 1.4616 and x0 +- 1e-5, 1e-10; tolerance
    4e-13 from the reflection formula's cancellation, 1313x at x = -22.75,
    where SciPy's error is 2.8e-13; tested) and `trigamma` (943 points, 1e-13;
    reference uses the reflection formula with `sinpi` for x < 0 because
    `mpmath.psi(1, x)` recurs one unit at a time from negative x; x =
    -999999999999999.5 left out because SciPy's `polygamma(1, x)` never
    returns there). SciPy passes every point of both.
  - `digamma` fails 7 points: **finding 017** (high), approved by Claude after
    the checklist. It recurses once per unit below 8, so it overflows the stack
    and aborts the process for large negative x (from about -1e5 in debug,
    -2e5 in release, on an 8 MiB stack; always below -2^53 and at -inf); and it
    loses all digits next to its zero (wrong sign at the nearest double, up to
    2.1e-3 away off the grid): plain cancellation, no expansion around x0 as
    in Boost and SciPy.
  - `trigamma` fails 31 points: **finding 018** (high), approved by Claude
    after the checklist. Finite at every negative-integer pole, and an error
    growing with |x| for negative x (4.8e-13 at -1000.25, 90x near -1e15):
    `sin(pi x)` of a rounded `pi x`. Confirmed: a `sinpi` removes every
    failure (0 of 5,111) and gives +inf at the poles.

- **Milestone 4, part 10: statrs Beta distribution** (statrs `Beta::cdf`,
  `sf`, `inverse_cdf`; [PR #22](https://github.com/traderhall757-arch/rust-sci-verify/pull/22), merged):
  - New datasets `beta_cdf` and `beta_sf` (1895 points each: the betainc
    a, b x 19 values of x from 0 to 1, including 1e-300 and the doubles
    nearest 1 - 1e-15 and 1, plus mean + k sd; tolerance 1e-12, betainc's;
    the largest prefactor logarithm at a normal result, ln Gamma(2e4) =
    1.78e5, would cost 2.0e-11 if rounded separately; tested). Reference: the
    continued fraction on the converging side with 1 - x exact, tested
    against `mpmath.betainc`, the series in both tails and closed forms.
    SciPy is wrong at a = b = 1/2 next to 0 (sf) and 1 (cdf). The quantile
    uses the existing `betaincinv` data.
  - `cdf` and `sf` fail 156 and 129 points: **finding 019** (high), approved
    by Claude after the checklist. `cdf` returns exactly 1 when x is within
    1e-15 of 1 (`beta_reg` sets its prefactor to 0 there; 9.1% off at
    a = 1e4, b = 0.1); both lose digits for max(a, b) >= 100 (prefactor from
    separately rounded logs; up to 5.0e-10 on the grid, 1.5e-9 off it).
    Confirmed: a bit-exact copy without the cut-off and with an exact
    prefactor fails no grid point.
  - `inverse_cdf` fails 479 of 1620: **finding 020** (high), approved by
    Claude after the checklist. It is special's `inv_inc_beta` (finding 015)
    with the same counts: 161 panics (p < 1.5e-162: NaN start, then
    `beta_reg` unwraps an out-of-range error, in release builds too),
    7 endless loops, 125 wrong small quantiles, 186 small errors.

- **Milestone 4, part 11: statrs F distribution** (statrs
  `FisherSnedecor::cdf`, `sf`, `inverse_cdf`;
  [PR #23](https://github.com/traderhall757-arch/rust-sci-verify/pull/23), merged):
  - New datasets `f_cdf` and `f_sf` (1999 points each: d1, d2 in {0.5, 1, 2,
    5, 10, 30, 100, 1000, 1e5}, every pair, x 21 values of x from 0 to 1e300
    plus 1 + k s around the bulk; tolerance 1e-12, betainc's; condition
    number at most 1648 and prefix logarithm 1.05e6, both tested) and
    `f_ppf` (1620 points: the same pairs x the betaincinv probabilities;
    p = 1 tagged `pole`; 2e-13 from the largest |log x| 690.8, tested).
    Reference: the continued fraction with z = d1 x / (d1 x + d2) and
    y = d2 / (d1 x + d2) formed separately; the quantile is (d2/d1) e^s with
    s = log(z / (1 - z)) from the betaincinv Newton iteration, factored out
    as `_mp_betaincinv_logit` (betaincinv data unchanged). Tested against
    `mpmath.betainc`, the series on the smaller tail, closed forms for
    d1 = d2 = 1, d1 = 2, d2 = 2, and F(ppf(p)) = p. SciPy fails 6 cdf,
    11 sf (up to 6.4e-11) and 101 ppf points (all with p <= 1e-100).
  - `cdf` and `sf` fail 194 and 440 points: **finding 021** (high), approved
    by Claude after the checklist. `sf` forms 1 - z from the rounded z, so it
    is 0 once d1 x / d2 exceeds about 2^53 (117 points) and loses digits as x
    grows (2.6% off at d1 = 1e5, d2 = 0.5, x = 1e10; 138% off the grid);
    `beta_reg`'s 1e-15 cut-off gives exactly 1 next to the ends (20 points
    each); finding 019's prefactor and the 140-iteration cap (which matters
    for shapes near 5e4) cost up to 1.3e-9. All four causes confirmed with a
    bit-exact switchable copy.
  - `inverse_cdf` fails 724 of 1620: **finding 022** (high), approved by
    Claude after the checklist. `inv_beta_reg`'s problems (finding 020):
    162 panics, 7 endless loops, 345 errors in z; plus 210 upper quantiles
    lost in `d2 / (d1 (1/z - 1))` even with the correctly rounded z (+inf,
    values capped near 2^52 d2 / d1, 0.32% off at the median for d1 = 1e5,
    d2 = 0.5).

- **Faster `reference-fresh`** ([PR #24](https://github.com/traderhall757-arch/rust-sci-verify/pull/24), merged):
  - `python -m rsv_reference` runs each dataset in its own fresh process, in
    parallel (driver in `rsv_reference/cli.py`): a full run takes 4 m 42 s
    locally instead of 14 m 46 s, with byte-identical data.
  - On pull requests, `scripts/select_datasets.py` regenerates only the
    datasets whose `Spec` reaches changed definitions in `functions.py`
    (syntax trees), whose data file changed, or that are new. It regenerates
    everything when the generator, environment, workflow or shared code
    changed, or when it cannot follow the code. It prints one line per
    dataset (demo: one Spec changed, 3 s). Guardrail 10 in `docs/guardrails.md`.
  - `main` (runs no longer cancelled) and a new weekly `reference-fresh` job
    (opens a `reference-stale` issue) still regenerate everything; the check
    also fails on untracked data files now.
- **Milestone 4, part 12: statrs Exponential and Erlang distributions**
  (statrs `Exp::cdf`, `sf`, `inverse_cdf`; `Erlang::cdf`, `sf`,
  `inverse_cdf`;
  [PR #25](https://github.com/traderhall757-arch/rust-sci-verify/pull/25), merged):
  - New datasets `expon_cdf` and `expon_sf` (325 points each: x = 10^-k for
    k = 1..20, 0..40 step 0.25, the underflow edge 700..750, decades from
    1e-300 to 1e300, negative x) and `expon_ppf` (138 points, the `norm_ppf`
    probabilities; p = 1 tagged `pole`). Tolerance 1e-13 for all three: the
    inputs are exact and expm1, exp and log1p of an exact argument are within
    an ulp. Reference: `-expm1(-x)`, `exp(-x)`, `-log1p(-p)`, tested against
    the gamma routes with shape 1, `mpmath.gammainc`, the Taylor series and
    `-log(1 - p)` at 1200 bits. SciPy passes every point (worst 2.2e-16).
  - `Exp::sf` is exact at every point; `Exp::inverse_cdf` passes (worst
    2.1e-16).
  - `Exp::cdf` fails 45 of 325: **finding 023** (high), approved by Claude
    after the checklist. It computes `1 - exp(-x)`, so it returns 0 for every
    x below 2^-54 = 5.6e-17 (1e-20 gives 0) and loses digits up to about
    x = 5e-4 (8.3e-8 at 1e-10, 5.7e-13 at 1e-5). Off the grid 9,805 of 20,000
    random x fail. Confirmed: `-expm1(-x)` passes every grid and random
    point.
  - `Erlang` calls `Gamma` for `cdf`, `sf` and `inverse_cdf`, so it is checked
    on the gamma datasets' whole-number shapes (409, 409 and 324 points; the
    report's grid line says so) and fails exactly where `Gamma` does:
    findings 012 and 013 now name it. `GAMMA_PPF` in `expected.rs` is split
    into `GAMMA_PPF_WHOLE` (111, shared with `Erlang`) and
    `GAMMA_PPF_FRACTIONAL` (59); the `Gamma` check lists both, unchanged.
  - Correction to the earlier plan: no statrs continuous distribution except
    `Dirac` uses the generic bisection `inverse_cdf` default (every other one,
    `Erlang` included, has its own), so there is nothing more to check there.
  - Reference generator fix found by this PR's CI: NumPy's AVX-512 `exp`
    differs from the C library's in the last bit, so `expon_sf`'s SciPy column
    depended on the CPU (9 points). `rsv_reference` now disables NumPy's
    `X86_V4` targets before NumPy is imported and refuses to run otherwise;
    regenerating every dataset changes only those 9 values.

- **Milestone 4, part 13: statrs InverseGamma distribution** (statrs
  `InverseGamma::cdf`, `sf`, `inverse_cdf`;
  [PR #26](https://github.com/traderhall757-arch/rust-sci-verify/pull/26), merged):
  - New datasets `invgamma_cdf` and `invgamma_sf` (531 points each: 15 shapes
    0.01 to 1e4 x x = 1/t for the gamma_cdf grid's t, plus 0, 5e-324 and
    1e-310; tolerance 1e-12, gamma's, since F(x) = Q(a, 1/x)) and
    `invgamma_ppf` (405 points: the gamma_ppf probabilities; p = 1 tagged
    `pole`, 3 `overflow`; 2e-13 from the largest |log x| 691.3, tested).
    Shapes stop at 1e4 so that rounding 1/x (which SciPy, Boost and statrs
    all do) costs at most 3.4e-13 (tested). a = 150 was added to the grid
    after the off-grid scan found the quantile's overflow window. The
    reference reuses the gamma routes at 1/x (60 extra bits); the gamma
    quantile was factored into `_mp_gamma_quantile` (either tail), with
    `gamma_ppf` byte-identical. Tested against `mpmath.gammainc` at 1/x and
    closed forms for a = 1 and 1/2. SciPy fails 1 cdf, 1 sf (flushed
    subnormals) and 11 ppf points (p = 5e-324).
  - `cdf` and `sf` **panic** for every x below 5.56e-309 (30 points each):
    **finding 024** (medium), approved by Claude after the checklist.
    `rate / x` overflows to inf and `gamma_ur` / `gamma_lr` unwrap an
    `XInvalid` error. Their other 8 and 10 failures (a >= 1000, worst
    2.5e-11, and two flushed subnormals) are finding 012's, now naming
    `InverseGamma`.
  - `inverse_cdf` fails 38 of 405: **finding 025** (high), approved by Claude
    after the checklist. For shapes from about 108 to 170 the density
    `x^(-a-1) e^(-1/x) / Gamma(a)` overflows to inf, the Newton step becomes
    zero and the search returns its first guess (median of
    `InverseGamma(150, 1)` 12% off; up to 46% off the grid); in the upper tail
    for small shapes the bracket comes from the rounded `cdf` and misses the
    root (up to 77% off the grid); and at subnormal p (up to 3.8%). Causes 1
    and 2 confirmed with a bit-exact copy of `newton_raphson_quantile`.
    `Chi::inverse_cdf` uses the same search.

- **Milestone 4, part 14: statrs one-sample t-test** (statrs
  `stats_tests::ttest_onesample`;
  [PR #27](https://github.com/traderhall757-arch/rust-sci-verify/pull/27), merged):
  - The generator now takes vector inputs (a sample) and several outputs
    (`docs/data-format.md`); every existing data file regenerates
    byte-identical.
  - New dataset `ttest_1samp` (324 points: deterministic samples of n = 2 to
    100 around 0 and 10, nine hypothesised means, three alternatives; outputs
    statistic and p-value; |t| from 0.017 to 3410, p down to 4.3e-253;
    tolerance 5e-12 from the rounding bound max(kappa_p, 1) (kappa_t + n)
    2^-53, largest 4.5e-12, tested). Reference: exact mean and variance, the
    `t_cdf` / `t_sf` routes; tested against exact rational arithmetic and the
    closed forms for 1 and 2 degrees of freedom. SciPy passes every point.
  - `ttest_onesample` fails 60 points: **finding 026** (high), approved by
    Claude after the checklist. Two-sided and "greater" p-values are
    `1 - cdf`: 0 below about 1e-16, digits lost below about 4e-5 (2.9% off
    at 1.8e-15). Off the grid also up to 2.3e-11 in t when the mean nearly
    equals mu (plain summation). Both confirmed: `sf` and a compensated sum
    leave 0 of 2,589 random samples failing.

- **Milestone 4, part 15: statrs chi-squared test and one-way ANOVA** (statrs
  `stats_tests::chisquare`, `stats_tests::f_oneway`;
  [PR #28](https://github.com/traderhall757-arch/rust-sci-verify/pull/28), merged):
  - New datasets `chisquare` (175 points: 2 to 30 categories, equal or
    sloping expected frequencies, totals 100 to 1e6, ddof 0 and 1) and
    `f_oneway` (48 points: 2 to 5 groups of 2 to 30 values around 0 and 100,
    group shifts 0.001 to 30); outputs statistic and p-value; tolerance 5e-12
    from the same rounding bound as `ttest_1samp` (largest 4.2e-12 in both,
    tested). Tested against exact rational arithmetic, closed forms,
    `mpmath.gammainc` / `mpmath.betainc` and, for two groups, the two-sample
    t-test. SciPy is off at 4 `f_oneway` points.
  - `chisquare` fails 29 points (14 of the equal-frequency points with
    `f_exp = None`): **finding 027** (high), approved by Claude after the
    checklist. The p-value is `1 - cdf`, as in finding 026.
  - `f_oneway` fails 34 of 48: **finding 028** (high), approved by Claude
    after the checklist. Textbook one-pass sums of squares lose digits when
    the values are large compared with their spread (F 1.7e-4 off for two
    groups of three values near 1e5), plus the `1 - cdf` p-value.
  - Reference generator fix found by this PR's CI: OpenBLAS picks its kernels
    by CPU, so SciPy's `f_oneway` differed in the last bit between a machine
    with AVX-512 (SkylakeX kernels) and the CI runners (Haswell). The
    generator now pins `OPENBLAS_CORETYPE=Haswell` before NumPy and SciPy are
    imported and refuses to run otherwise; regenerating every dataset changes
    only 15 SciPy values of `f_oneway`.

- **Milestone 4, part 16: statrs skewness test** (statrs
  `stats_tests::skewtest`;
  [PR #29](https://github.com/traderhall757-arch/rust-sci-verify/pull/29), merged):
  - New dataset `skewtest` (210 points: samples of 8 to 100 values with
    quadratic skew, an outlier, or exact symmetry, around 0 and 10; three
    alternatives; outputs z and p-value; tolerance 2e-12 from the rounding
    bound, largest 1.6e-12, tested). Reference: D'Agostino's transform from
    the exact sample skewness; tested against exact rational arithmetic and
    the published SciPy example.
  - `skewtest` fails 149 points. **Finding 029** (high), approved by Claude
    after the checklist, covers 49: it replaces a skewness of exactly 0 by
    1 (as SciPy does), so symmetric data such as 1..8 get z = 1.01 and a
    two-sided p of 0.31 instead of 1; `1 - cdf` p-values; plain summation.
    The other 100 fail only through `Normal::cdf`'s `erfc` (finding 009,
    now naming `skewtest`).

- **Milestone 4, part 17: statrs ln_gamma and digamma** (follow-up; statrs
  `function::gamma::ln_gamma`, `digamma`, on the existing `gammaln` and
  `digamma` data;
  [PR #30](https://github.com/traderhall757-arch/rust-sci-verify/pull/30), merged):
  - `ln_gamma` is documented as the logarithm of Gamma, so it runs on the
    877 `gammaln` points where Gamma(x) > 0 or x is a pole (its NaN where
    Gamma < 0 is its documented answer). It fails 32: **finding 030**
    (high), approved by Claude after the checklist. Finite values at 19
    negative-integer poles (`ln_gamma(-2)` = 36.4), and all digits lost next
    to 1 and 2 (`ln_gamma(1)` = -4.4e-16, `ln_gamma(1 - 1e-15)` = 0): a sum of
    rounded logarithms that cancel, and a sine of a rounded pi x.
  - `digamma` fails 9 of 949: **finding 031** (high), approved by Claude after
    the checklist. Next to its zero (14x off at the nearest double) and at
    large negative x (1.5% off at -1e15 - 0.5, from tan of a rounded pi x;
    confirmed with an exactly reduced argument). No stack overflow, unlike
    special's (finding 017).

- **Milestone 4, part 18: statrs gamma_lr and beta_reg** (follow-up; statrs
  `function::gamma::gamma_lr` on `gammainc` and `gammainc_large`,
  `function::beta::beta_reg` on the three betainc datasets;
  [PR #31](https://github.com/traderhall757-arch/rust-sci-verify/pull/31), merged):
  - `gamma_lr` panics at x = 0 for every shape (12 grid points; `gamma_ur`
    too): **finding 032** (medium), approved by Claude after the checklist.
    Its other 15 failures are finding 012's (now naming it): up to 4.1e-6 at
    a = 1e9, the spacing of doubles near a ln a, and one flushed subnormal.
  - `beta_reg` fails 195 of 969 points: finding 019's prefactor (now naming
    it), plus statrs's 140-step cap at a = b = 1e5; a bit-exact copy with
    the exact prefactor and no cap passes all 969.

- **Milestone 4, part 19: puruspe ln_gamma and gammp for large shapes**
  (follow-up; `puruspe::ln_gamma` on `gammaln`, `puruspe::gammp` on
  `gammainc_large`;
  [PR #32](https://github.com/traderhall757-arch/rust-sci-verify/pull/32), merged):
  - `ln_gamma` runs on the 877 `gammaln` points where Gamma(x) > 0 or x is
    a pole (as for statrs) and fails 666: **finding 033** (high), approved
    by Claude after the checklist. No reflection formula (NaN below -4.5,
    up to 1.8e-4 off on (-4, -3) and (-2, -1)); +inf for x below 2^-54
    (x - 1 rounds to -1); `ln_gamma(1)` = -5.2e-12, so all digits are lost
    next to 1 and 2; and an absolute error up to 1.3e-10 everywhere (its
    7-term Lanczos set), above the tolerance up to x = 156.
  - `gammp` fails 25 of 51 `gammainc_large` points: finding 002, extended
    and now root-caused. NaN at all 19 grid points with 0.2a <= x <= a - 1:
    puruspe's log-sum-exp rewrite of Numerical Recipes' `gammpapprox` takes
    the logarithm of the negative `xu - x` (the same cause as the original
    NaN at a >= 100; `gammq` too). Up to 1.5e-6 near the mean (1.3e-5 off
    the grid at a = 1e10): the a ln a cancellation in the scale and the
    18-point quadrature, of equal size (shown with a copy of the function).

- **Milestone 4, part 20: puruspe gamma, gammq, and gammp on the
  distribution data** (follow-up; `puruspe::gamma` on `gamma`,
  `puruspe::gammq` on `gamma_sf`, `puruspe::gammp` on `gamma_cdf`;
  [PR #33](https://github.com/traderhall757-arch/rust-sci-verify/pull/33), merged):
  - `gamma` (the Fukushima port) passes all 901 points (worst 1.5e-15).
  - `gammq` fails 223 of 538. **Finding 034** (high), approved by Claude after
    the checklist, covers 180: it returns 0 for every x > a + 20 (sqrt(a) + 1)
    (63 points; `gammq(1, 50)` = 0 instead of e^-50 = 1.9e-22), 1 for every
    x < a 2^-52 (6 points; `gammq(0.01, 1e-20)` = 1 instead of 0.365), and
    loses digits through `ln_gamma` (finding 033) in its prefactor (111
    points, up to 8.1e-9 where Q = 1 - P is small). The other 43 are the
    shared quadrature (finding 002: 23 NaN, 20 up to 1.2e-9).
  - `gammp` on `gamma_cdf` fails 175 of 538: 20 from the same small-x
    shortcut (returns 0; `gammp(0.01, 1e-20)` = 0 instead of 0.635; finding
    034), which the `gammainc` grid (x >= 0.001) never reached; 155 are
    finding 002's (extended).
  - Finding 002's accuracy loss (group 2) and flushed subnormal (group 3)
    are now confirmed: a bit-exact switchable copy with libm's `lgamma` and
    without the `log_result < -700` shortcut leaves only the quadrature
    region failing on `gammainc`, `gamma_cdf` and `gamma_sf`.

- **Milestone 4, part 21: statrs Chi distribution** (statrs `Chi::cdf`,
  `sf`, `inverse_cdf`;
  [PR #35](https://github.com/traderhall757-arch/rust-sci-verify/pull/35), merged):
  - New datasets `chi_cdf` and `chi_sf` (514 points each: 13 degrees of
    freedom k from 1 to 1e4, including 160 and 161 around statrs's switch of
    density formula, x x = sqrt(2 t) for the gamma_cdf grid's t at shape
    k/2, plus x where x * x / 2 underflows or x * x overflows; tolerance
    1e-12, gamma's, since F(x) = P(k/2, x^2/2); rounding x * x costs at most
    2.4e-13, tested) and `chi_ppf` (351 points: the gamma_ppf probabilities;
    p = 1 tagged `pole`; 2e-13 from the largest |log(x / sqrt 2)| 690.9,
    tested). The reference reuses the gamma routes at x^2/2 (60 extra bits)
    and the gamma quantile; tested against `mpmath.gammainc`, closed forms
    for k = 1, 2 and 3, and the finite sum at k = 1e4. SciPy fails 5 cdf, 4
    sf and 4 ppf points.
  - `cdf` and `sf` **panic** wherever `x * x / 2` rounds to 0 (x below
    2.72e-162) or `x * x` overflows (from 1.34e154), at every k (52 points
    each), and for k = 1 `cdf` is 0 for x up to 7.2e-162 (`gamma_lr` treats
    the five smallest subnormals as 0) and up to 3.7e-4 off up to about
    x = 2e-156 (a subnormal square): **finding 035** (high), approved by
    Claude after the checklist. Their other 6 and 7 failures (k = 1e4, worst
    8.6e-12, and flushed subnormals) are finding 012's, now naming `Chi`.
  - `inverse_cdf` fails 14 of 351: **finding 036** (high), approved by
    Claude after the checklist. For k = 1 every p below about 6e-162 gives
    6.67e-162 or 7.20e-162 (`inverse_cdf(1e-300)` instead of 1.25e-300),
    because the cdf is 0 there; at p = 5e-324 every k is wrong (a factor of
    3e7 at k = 2, 10% at k = 300), because the cdf flushes values below
    1e-308 to 0 (finding 012) and then finding 025's cause 3. Unlike
    `InverseGamma`, no density overflow and no bracket miss (5,000 + 3,000
    point off-grid scan, seed 20261016). Confirmed with a bit-exact copy of
    the search.

- **Milestone 4, part 22: statrs Weibull distribution** (statrs
  `Weibull::cdf`, `sf`, `inverse_cdf`;
  [PR #36](https://github.com/traderhall757-arch/rust-sci-verify/pull/36), merged):
  - New datasets `weibull_cdf` and `weibull_sf` (1139 points each: shapes
    0.1 to 100 x scales 2^-500, 2^-10, 1, 2^10, 2^500, x = u * scale with
    u^k from 1e-300 to 1e300, plus x = 0 and -scale) and `weibull_ppf`
    (1215 points, the gamma_ppf probabilities; p = 1 tagged `pole`).
    Tolerance 1e-13. Scales are powers of two so x / scale is exact (the
    sf's condition number in x is k t, up to about 70,000, so otherwise no
    method could meet 1e-13); this leaves non-power-of-two scales
    unchecked. Reference `-expm1(-t)`, `exp(-t)` with t = (x/scale)^k, and
    scale (-log1p(-p))^(1/k), tested against the exponential (shape 1),
    quadrature of the density and cdf(ppf(p)) = p. SciPy passes every point.
  - All three fail where statrs's stored `scale.powf(-shape)` (or
    `x.powf(shape)`, or `-ln_1p(-p) / scale.powf(-shape)`) overflows,
    underflows or is subnormal although the answer is moderate: NaN, 0, 1
    or inf. `cdf` 177 and `sf` 180 of 1139: **finding 037** (high);
    `inverse_cdf` 214 of 1215: **finding 038** (high); both approved by
    Claude after the checklist. Ordinary parameters are hit too (shape 100,
    scale 1/1024: `cdf` wrong for every x below 7.8e-4). Random scans (seed
    20261017): every failure but one sf rounding case has an out-of-range
    intermediate; computing t = (x / scale)^shape first (and multiplying the
    quantile by scale last) fixes every grid point.
  - `pdf` and `ln_pdf` use the same stored power; not checked.

- **Milestone 4, part 23: statrs log-normal distribution** (statrs
  `LogNormal::cdf`, `sf`, `inverse_cdf`;
  [PR #37](https://github.com/traderhall757-arch/rust-sci-verify/pull/37), merged):
  - New datasets `lognorm_cdf` and `lognorm_sf` (652 points each: (mu,
    sigma) in {(0, 0.25), (-3, 0.5), (0, 1), (2, 4)} x x = exp(mu + sigma w)
    for w = -40..40 step 0.5, plus x = 0 and -1; tolerance 8e-13, because
    rounding ln x and w costs up to 3.4e-13 in the tails) and `lognorm_ppf`
    (552 points, the norm_ppf probabilities; p = 1 tagged `pole`; tolerance
    1e-13). Reference Phi((ln x - mu) / sigma) and exp(mu + sigma
    Phi^-1(p)), tested against direct quadrature of the density in both
    tails (out to |w| = 37.5, no erfc involved), erfinv and cdf(ppf(p)) = p;
    a comparison with `mpmath.ncdf` checks only the wiring (it uses the
    same erfc for w < 0). SciPy passes
    every point except the subnormal results at w = -38 and 38 (flushed to
    0, as for the normal).
  - `cdf` and `sf` fail 308 of 652 each (about 10 digits, up to 8.6e-11,
    for 1 <= |w| <= 37.5 on the error-carrying side): **finding 039** (low),
    approved by Claude after the checklist; root cause statrs's `erfc`
    (finding 004), exactly as for `Normal` (finding 009). A copy with
    finding 004's fix passes every grid point and 40,000 random points
    (seed 20261018). `inverse_cdf` passes all points (worst 3.3e-14).
  - `scripts/tests/test_select_datasets.py`: the datasets that depend on
    `_mp_erfcx` now include the three log-normal ones (a guardrail test's
    expected list, listed under "Checks changed").

- **Milestone 4, part 24: statrs Pareto distribution** (statrs
  `Pareto::cdf`, `sf`, `inverse_cdf`;
  [PR #38](https://github.com/traderhall757-arch/rust-sci-verify/pull/38), merged):
  - New datasets `pareto_cdf` and `pareto_sf` (1062 points each: shapes
    2^-5, 0.1 to 100 x scales 2^-500, 2^-10, 1, 2^10, 2^500, x = u x_m with
    u from 1 + 2^-52 to 1e300 plus 1.1 2^1030, 1.1 2^1060, 2^1100 and
    2^1500 where x is finite, plus four points below the support) and
    `pareto_ppf` (1215 points, the gamma_ppf probabilities; p = 1 tagged
    `pole`, 10 results tagged `overflow`). Tolerance 1e-13. Scales are
    powers of two, so x / x_m is exact; shape 2^-5 makes 1/alpha = 32
    exact. Reference -expm1(-t), exp(-t) with t = alpha log(x / x_m), and
    x_m exp(-log1p(-p) / alpha), tested against the binomial series near
    the scale, quadrature of the density, the exact far tail and
    cdf(ppf(p)) = p (no log or expm1 in these); the exponential test checks
    the log only. SciPy is also wrong at 119 cdf, 18 sf and 1 ppf points
    (the same two causes as statrs).
  - `cdf` fails 121 and `sf` 9 of 1062: **finding 040** (high), approved by
    Claude after the checklist. Cause 1 (cdf): 1 - (x_m/x)^alpha cancels
    just above the scale (0 instead of 6.9e-18; about 10 digits at shape 1,
    scale 1, x = 1 + 1e-10). Cause 2: x_m / x is subnormal or 0 although the
    answer is moderate for shapes below 1 (sf 0 instead of 4.5e-11).
  - `inverse_cdf` fails 1 of 1215: **finding 041** (high, narrow range),
    approved by Claude after the checklist: (1 - p)^(-1/alpha) overflows
    for shapes below about 0.05 although x_m times it is finite.
  - Random scan (seed 20261019, 40,000 points, scales not powers of two):
    cdf 7,937, sf 989, inverse_cdf 44 failures, all from these causes; a
    copy with remedies (log1p of the exact x - x_m, scaling x_m by 2^1000,
    half the power applied twice) passes every grid and random point.
  - `pdf` and `ln_pdf` not checked.

- **Milestone 4, part 25: statrs Gumbel distribution** (statrs
  `Gumbel::cdf`, `sf`, `inverse_cdf`;
  [PR #39](https://github.com/traderhall757-arch/rust-sci-verify/pull/39), merged):
  - New datasets `gumbel_cdf` and `gumbel_sf` (1092 points each: (mu,
    beta) in (0, 1), (0, 2^-500), (0, 2^500), (-10, 1/4), (100, 8),
    (-2^1023, 2^1023); x = mu + beta z for z = -7..40 step 1/4 plus far
    tails and +-2^-30, kept where x is exact, so z is exact) and
    `gumbel_ppf` (204 points, the gamma_ppf probabilities plus the doubles
    nearest 1/e and 1/e +- 10^-k, k = 5, 10, 15; p = 0 and 1 tagged
    `pole`, 15 results tagged `overflow`). Tolerance 1e-13. Reference
    exp(-w), -expm1(-w) with w = exp(-z), and mu - beta log(-log p),
    tested against quadrature of the density, the upper-tail series and
    the root of cdf(x) = p at 1200 bits. SciPy is wrong at the same points
    as statrs.
  - `cdf` and `sf` fail 4 of 1092 each: **finding 042** (high, narrow
    range), approved by Claude after the checklist: x - mu overflows when
    both are near the largest double (cdf 1 instead of 0.87). Off the grid
    `cdf` also loses up to 3.2e-13 in the far lower tail from rounding z
    once (noted in the finding, not counted as a bug).
  - `inverse_cdf` fails 20 of 204: **finding 043** (high), approved by
    Claude after the checklist: log(-log p) cancels next to p = 1/e (0
    instead of 3.4e-17; 3.8e-13 off at scale 10, p = 0.3679), and beta
    log(-log p) overflows although the quantile is finite.
  - Random scans (seeds 20261020 and 20261021): all failures from these
    causes, except 2 quantile points where mu - beta L itself cancels
    (any method that rounds L loses there); a copy with remedies passes
    every grid and random point otherwise.
  - `pdf` and `ln_pdf` not checked.

- **Milestone 4, part 26: statrs Laplace distribution** (statrs
  `Laplace::cdf`, `sf`, `inverse_cdf`;
  [PR #40](https://github.com/traderhall757-arch/rust-sci-verify/pull/40), merged):
  - New datasets `laplace_cdf` and `laplace_sf` (1828 points each: the
    Gumbel (mu, b) pairs x x = mu + b z for z = -40..40 step 1/4, both far
    tails to +-1e300 and +-2^-30, +-2^-60, kept where x is exact) and
    `laplace_ppf` (174 points, the gamma_ppf probabilities plus the doubles
    next to 1/2; p = 0 and 1 tagged `pole`, 16 `overflow`). Tolerance
    1e-13. Reference e^z / 2 and 1 - e^-z / 2, and mu + b log(2p) /
    mu - b log(2 - 2p), tested against quadrature of the density, half the
    exponential tail and the root of cdf(x) = p at 1200 bits. SciPy is
    wrong at the same overflow points as statrs.
  - `cdf` and `sf` fail 4 of 1828 each: **finding 044** (high, narrow
    range), approved by Claude after the checklist: x - mu overflows
    (finding 042's cause), giving exactly 1 and 0 (or 0 and 1).
  - `inverse_cdf` fails 13 of 174: **finding 045** (high), approved by
    Claude after the checklist: it **panics** at p = 0 and p = 1 for every
    parameter (12 points), and b log(2 - 2p) overflows although the
    quantile is finite (1 point).
  - Random scans (seeds 20261022 to 20261024): all failures from these
    causes; a copy with remedies passes every grid and random point.
  - `pdf` and `ln_pdf` not checked.

- **Milestone 4, part 27: statrs Cauchy distribution** (statrs
  `Cauchy::cdf`, `sf`, `inverse_cdf`;
  [PR #41](https://github.com/traderhall757-arch/rust-sci-verify/pull/41), merged):
  - New datasets `cauchy_cdf` and `cauchy_sf` (1762 points each: the
    Gumbel (mu, b) pairs x x = mu + b z for z = -40..40 step 1/4, both
    heavy tails to +-1e308 (subnormal results) and +-2^-30, +-2^-60, kept
    where x is exact) and `cauchy_ppf` (186 points, the gamma_ppf
    probabilities plus the doubles next to 1/2 and p = 1/8, 7/8; p = 0 and
    1 tagged `pole`, 25 `overflow`). Tolerance 1e-13. Reference acot(-z) /
    pi and 1 - acot(z) / pi, and mu - b cospi(p) / sinpi(p), tested
    against quadrature of the density, the acot series and the root of
    cdf(x) = p at 1200 bits. SciPy is right in the tails, wrong at the
    overflow points.
  - `cdf` and `sf` fail 56 of 1762 each: **finding 046** (high), approved
    by Claude after the checklist: atan(z) / pi + 0.5 cancels in the tail
    (0 instead of 3.2e-17 at z = -1e16; 2e-7 off at z = -1e10), and
    x - mu overflows (finding 044's cause).
  - `inverse_cdf` fails 62 of 186: **finding 047** (high), approved by
    Claude after the checklist: finite at p = 0 and 1 (-1.6e16 b), rounding
    of pi (p - 1/2) next to 0 and 1 (every p below 2^-54 gives the p = 0
    value; 5e-7 off at p = 1e-10), tan(pi / 4) != 1 (nonzero where the
    quantile is 0), and b tan overflowing although the quantile is finite.
  - Random scans (seeds 20261025 to 20261027): all failures from these
    causes; a copy with remedies passes every grid and random point.
  - `pdf` and `ln_pdf` not checked.

- **Milestone 4, part 28: statrs Levy distribution** (statrs
  `Levy::cdf`, `sf`, `inverse_cdf`;
  [PR #42](https://github.com/traderhall757-arch/rust-sci-verify/pull/42), merged):
  - New datasets `levy_cdf` and `levy_sf` (1776 points each: the Gumbel
    (mu, c) pairs x x = mu + c t for t from -1 to 160/4 in steps of 1/65536
    to 1/4 next to the support, decades to 1e308, and 3 2^1030 to 2^1523,
    which only scale 2^-500 reaches; kept where x is exact) and `levy_ppf`
    (180 points, the gamma_ppf probabilities plus the doubles next to 1/2
    and 1; p = 1 tagged `pole`, 10 `overflow`; p = 0 gives mu). Tolerance
    1e-13. Reference erfc(s), erf(s) with s = sqrt(c / (2 (x - mu))), and
    mu + c / Phi^-1(p / 2)^2, tested against quadrature of the density, the
    erf series in the far tail and the root of cdf(x) = p at 1200 bits.
  - `cdf` fails 873 and `sf` 381 of 1776: **finding 048** (high), approved
    by Claude after the checklist: finding 004's erfc constants (erf is
    1 - erfc inside statrs, so both carry about 1e-10), erfc of a rounded
    square root next to the support (up to 1.5e-13 even with 004 fixed),
    x - mu overflow (exactly 1 and 0), and c / (2 (x - mu)) underflowing in
    the far upper tail (`sf` 0 instead of 4.7e-230).
  - `inverse_cdf` fails 5 of 180: **finding 049** (high, narrow range),
    approved by Claude after the checklist: c / (2 y^2) overflows although
    the quantile is finite (location -2^1023, p near 1/2).
  - Random scans (seeds 20261028 to 20261032): all failures from these
    causes, apart from one quantile point where mu and the tail term cancel
    (1815x, unavoidable); a copy with remedies passes every grid point and
    every other random point.
  - `pdf` and `ln_pdf` not checked.

- **Milestone 4, part 29: statrs Mann-Whitney U test** (statrs
  `stats_tests::mannwhitneyu`;
  [PR #43](https://github.com/traderhall757-arch/rust-sci-verify/pull/43), merged):
  - New dataset `mannwhitneyu` (690 points: two samples from the distinct
    integers (7919 i) mod 1000, x shifted by d; sizes 3 + 6 to 20 + 8 with
    all four methods (automatic, exact, asymptotic with and without
    continuity correction), 9 + 9 to 25 + 10 asymptotic only, and tied
    samples; three alternatives; outputs U1 and p-value; tolerance 1e-13 from
    the rounding bound, largest 8.6e-14, tested). Reference: exact average
    ranks in integer arithmetic, the exact null distribution from the Mann
    and Whitney recurrence, the tie-corrected normal tail from erfc; tested
    against enumeration of every split, pair counting, the exact variance
    with ties, `mpmath.ncdf` and SciPy's published example. SciPy passes
    every point (worst 2.1e-15).
  - `mannwhitneyu` fails 448 points. **Finding 050** (high), approved by
    Claude after the checklist, covers 180: the exact p-value is
    1 - P(U >= u) whenever x is not longer than y (SciPy's own example gives
    0.944 instead of 0.0556; every such exact point fails), and the
    asymptotic p-value is `1 - cdf` (0 below about 1e-16). The other 268
    fail only through `Normal::cdf`'s `erfc` (finding 009, now naming
    `mannwhitneyu`).
  - Random scan (seed 20261029, 3,000 points): all failures from these
    causes; a copy with both remedies and an accurate erfc passes every
    grid and random point. Also measured: the exact method enumerates every
    split (2.5 s for samples of 8 and 40), and `Automatic` picks it whenever
    one sample has at most 8 values (in the finding's notes, no listed
    points).

- **Milestone 4, part 30: statrs Fisher's exact test** (statrs
  `stats_tests::fisher::fishers_exact_with_odds_ratio`;
  [PR #44](https://github.com/traderhall757-arch/rust-sci-verify/pull/44), merged):
  - New dataset `fisher_exact` (411 points: published tables (SciPy's
    documentation, Fisher's tea tasting, statrs's fisher test cases), the first 12
    tables with another table at most 1e-4 more probable, and tables by their
    margins from 8 to 1900 values with the first cell from end to end of its
    support; three alternatives; outputs odds ratio and p-value; tolerance
    1e-13 from the rounding bound, largest 7.8e-14, tested). Reference:
    exact table counts C(n1, x) C(n2, n - x) in integer arithmetic, tested
    against enumeration of labellings, `mpmath.binomial` at 60 digits and
    the published values. SciPy passes every point (worst 2.0e-15). Tables
    with b c = 0 (odds ratio +inf) are left out: the generator would tag
    them `overflow`, and a finite p-value there is judged crate wrong
    (statrs's p-value is wrong on them too; in finding 051's notes).
  - `fishers_exact_with_odds_ratio` fails 165 points. **Finding 051** (high),
    approved by Claude after the checklist, covers 51 two-sided points: the
    upper tail is `1 - cdf` (2.2e-14 instead of 2.2e-59, and a negative
    p-value for a table of statrs's own tests), and a 1e-4 tie window adds
    more probable tables (0.529 instead of 0.368). **Finding 052** (high),
    approved likewise, covers the other 114: `Hypergeometric::pmf` is NaN
    once a binomial coefficient exceeds the largest double (populations
    from about 1030; a two-sided p of 0.493 instead of 0.051), and
    `Hypergeometric::cdf` sums exp of large log terms (one-sided p-values up
    to 1.5e-12 off on the grid, 6.1e-12 off it).
  - Random scan (seed 20261030, 3,000 tables) and a search for near ties:
    all failures from these causes; a bit-exact copy with all four remedies
    passes every grid, random and near-tie point.

- **Milestone 4, part 31: statrs one-sample Kolmogorov-Smirnov test** (statrs
  `stats_tests::ks_test::ks_onesample` against `Uniform(0, 1)`;
  [PR #45](https://github.com/traderhall757-arch/rust-sci-verify/pull/45), merged):
  - New dataset `ks_1samp` (544 points: scrambled samples z_i =
    ((7919 (i + 1)) mod 2048 + 1/2) / 2048 pushed towards 0 or 1 by c in
    {31/32 .. 1/2}, and evenly spread samples (2 i - 1) / (2 n), for n from
    1 to 150 with all five methods; n = 1100 without the exact one; outputs
    statistic and p-value; tolerance 5e-12 from the rounding bound,
    largest 3.3e-12, tested). Reference: the statistic exactly in rationals;
    the one-sided p-value from the Birnbaum and Tingey sum in exact integer
    arithmetic; two-sided exact from 2 P(D+ >= d) for d >= 1/2 and the
    Marsaglia, Tsang and Wang matrix below; the limit distribution from its
    two series. Tested (to 38 digits) against Steck's determinant in exact rationals (n up
    to 150), the published K(10, 0.274), closed forms, `mpmath.jtheta` and
    counting. SciPy is wrong at 4 exact points at n = 150 (up to 2.6e-6; it
    approximates above n = 140).
  - `ks_onesample` fails 30 points: **finding 053** (high), approved by
    Claude after the checklist. Four causes: the exact two-sided p-value is
    `1 - cdf` (2.7e-15 instead of 6.6e-27); `n.powi(n)` overflows from
    n = 144, so the exact method gives NaN (or a spurious 1) for 144 to 169
    values; the
    one-sided and approximate p-values overflow C(n, j) from n = 1030 (1 or
    NaN; 1 instead of 0.117) and underflow factors from about n = 700 when
    p < 1e-89; the asymptotic series stops at a term below 1e-10 (up to
    6.8e-11 off where p = 1).
  - Random scans (seeds 20261031, 3,000 samples, and 20261101, 1,000 nearly
    even samples): all failures from these causes; a copy with remedies for
    causes 2 to 4 passes everything except cause 1 below d = 1/2 (no simple
    remedy; absolute error at most 4e-15) and 2 points below 1e-303.
  - Guardrail change: `scripts/select_datasets.py` allows the standard
    library's `fractions` import in `functions.py` (exact rationals, no state).

- **Milestone 4, part 32: statrs two-sample Kolmogorov-Smirnov test** (statrs
  `stats_tests::ks_test::ks_twosample`; PR open, see "Next"):
  - New dataset `ks_2samp` (363 points: the Mann-Whitney samples with
    half-integer shifts and an interleaved pair, sizes (1, 1) to (100, 100)
    with all four methods and (101, 100) without the exact one; tied samples
    and y = x repeated (D = 0) with methods 0, 1 and 3; outputs statistic
    and p-value; tolerance 5e-12 from the rounding bound, largest 4.6e-12,
    tested). Reference: the statistic exactly in rationals; the exact
    p-value by counting lattice paths in exact integers (tested against
    enumeration and the Gnedenko and Korolyuk sum); Hodges's eq. 5.3 and
    the limit distribution as published formulas. SciPy agrees everywhere
    except a statistic of 1.1e-16 instead of 0 at 3 tied points.
  - `ks_twosample` fails 30 points: **finding 054** (high), approved by
    Claude after the checklist. Four causes: the empirical cdfs are summed
    1/n at a time, so a statistic that should be 0 is 1e-16 to 3e-16;
    the two-sided asymptotic method **never returns** when the statistic is
    0 or that rounding noise (samples with equal empirical cdfs, which need
    ties); the exact p-value is `1 - cdf` (-6.0e-15 instead of 2.2e-59);
    the asymptotic series stops at 1e-10 (finding 053's cause 4).
  - Random scan (seed 20261102, 3,000 pairs): all 330 failures from these
    causes (102 never return); a copy with remedies for all four passes
    everything.
  - The check runs first in the statrs test, in child processes
    (`run_isolated`), so the 4 endless calls are killed instead of left
    spinning; this adds about 20 s (4 x the 5 s limit) to the test.

## Next

1. **Owner: review and merge the open PR** (Milestone 4 part 32, statrs
   `ks_twosample`, finding 054).
2. **Milestone 4, continued** (next session; one family per PR):
   - **Remaining statrs hypothesis test** (suggested next):
     `anderson_darling` (its p-value is an empirical approximation by
     design: check the statistic, and the p-value only against the same
     published formula).
   - puruspe `invgammp` (calls `gammp`, so finding 034's small-x shortcut
     and finding 002's NaN may break it) on `gamma_ppf`; puruspe `invbetai`
     on `betaincinv`; puruspe `erf`/`erfc`/`inverfc` on `erfc`/`norm_ppf`.
   - Not yet listed anywhere: statrs `Triangular` and `Uniform` (closed form,
     continuous) and the discrete distributions (`Poisson`, `Binomial`,
     `Geometric`, `Hypergeometric`, ...); the owner decides whether they belong
     in Milestone 4. `Hypergeometric` already has finding 052 (seen through
     Fisher's test).
   - Optional: the K result of `Inu_Knu` on the `kv` grid (finding 008's
     code path).
   - Keep expected lists small: vector inputs are listed value by value, so
     samples stay at n <= 100 (see PR #29).
3. Nothing is submitted upstream until the whole project is done and tested
   (owner, 2026-10-14). The puruspe fork (`traderhall757-arch/puruspe`,
   branch `fix-betai`) still holds the old commit; replace it with the
   revised patch only when submission is approved.

## Decisions recorded (owner, 2026-10-14)

- Findings 001 to 004 approved; finding 005 approved (same day, after PR #12).
  Finding 006 and the J part of 005 approved by Claude (checklist); findings
  007 to 054 likewise (and finding 002's two extensions, 2026-10-15).
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
- `t_ppf` takes about 3 minutes to generate (bisection plus Newton on the
  continued fraction, at 50 and 100 digits) and `t_cdf` + `t_sf` about
  7 s; `gamma_ppf` about 5.5 minutes (2 of them at a = 1e6) and
  `gamma_cdf` + `gamma_sf` about 10 s; `betaincinv` about 5.5 minutes;
  `gammainc_large` about 1.5 minutes, `gammaln` and `gammasgn` about 2 s,
  `digamma` and `trigamma` about 3 s together; `beta_cdf` and `beta_sf`
  about 17 s together; `f_cdf` and `f_sf` about 23 s together, `f_ppf`
  about 6 minutes; the three `expon` datasets about 3 s; `chi_ppf` about
  15 s and `chi_cdf` + `chi_sf` about 1 s; the three `weibull` datasets under 1 s; the three `gumbel` datasets about 1 s; the three `laplace` datasets about 1 s; the three `levy` datasets about 1 s; `mannwhitneyu` under 1 s; `fisher_exact` about 2 s; `ks_1samp` about 8 s; `ks_2samp` under 1 s;
  a full regeneration now takes about 26 minutes.
- `yv` and `yv_large` take about 10 s together; `jv` and `jv_large` about 1 s; the three `iv` datasets about 1 s; the three `norm` datasets about 2 s. The `yv` grids leave out
  x = 0 (a pole; `Jnu_Ynu` documents a panic there).
- `betainc_large` takes about a minute to generate (mpmath series at a, b = 1e5);
  `betainc_switch` about 30 s. Its grid leaves out a = 1e5 with x = 1 - 1e-4,
  where the reference series needs millions of terms.
- Random spot-checks of `betai` are done with scratch scripts (not committed,
  rule 6); the method and seed are recorded in finding 001.
- Upstream puruspe `main` fails its own `cargo fmt --check` and `cargo clippy`
  on Rust 1.97, independently of our patch (`fixes/README.md`).
- The reference generator pins NumPy below AVX-512 and OpenBLAS to its
  Haswell kernels (`rsv_reference/__init__.py`), so SciPy columns do not
  depend on the CPU; it refuses to run otherwise. Run the whole reference
  `pytest` (not one file): a single test file that imports NumPy first
  defeats the pins (`test_generate.py` alone fails its two pin tests).
- Hypothesis-test tolerances follow one rule: the largest bound
  max(kappa_p, 1) (kappa + n) 2^-53 a careful method can reach on the grid,
  rounded up (tested per dataset).
- The weekly job has been dry-run locally; its first real run is the first
  Monday after merge, or run it by hand from the Actions tab.
