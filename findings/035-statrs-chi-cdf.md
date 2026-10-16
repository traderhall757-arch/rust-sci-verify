# 035: statrs Chi cdf and sf panic where x * x under- or overflows, and are wrong for k = 1 below x = 2e-156

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Chi`, methods
  `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: `Chi::new(1).cdf(4e-162)` is
  0 instead of 3.19e-162, and `cdf(8.7e-161)` is 3.7e-4 off; plus panics,
  which end the calling thread, for every k where `x * x` under- or
  overflows)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 54 of 514 points of `reference/data/chi_cdf.jsonl`
  and 52 of 514 of `reference/data/chi_sf.jsonl` (tolerance 1e-12). The
  exact lists are `CHI_SQUARE_RANGE` (both methods) and
  `CHI_CDF_SUBNORMAL_SQUARE` (`cdf` only) in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Three groups, by
  cause (section "Root cause"):
  1. **Panics** (52 points of each dataset: x = 5e-324, 1e-170, 1e155 and
     1e300 at every one of the 13 k): every x below 2.72e-162 (where
     `x * x / 2` rounds to 0) and every x from 1.34e154 up (where `x * x` is
     +inf), at every k. The true values there are 0 or 1, except `cdf` for
     k = 1, which is about 0.8 x (7.98e-171 at x = 1e-170).
  2. **k = 1, `x * x / 2` among the five smallest subnormals** (1 point,
     `cdf` at x = 4e-162): every x from 2.72e-162 to 7.20e-162 gives 0
     instead of about 0.8 x.
  3. **k = 1, `x * x / 2` subnormal** (1 point, `cdf` at x = 1e-160,
     5.6e-6 off): the subnormal square keeps only a few bits, so `cdf` is
     off by up to 3.7e-4 for x up to about 2e-156 (scan).

The other `Chi::cdf` and `sf` failures (6 and 7 points: k = 1e4, worst
8.6e-12, and subnormal results returned as 0) are finding 012's errors in
`gamma_lr` and `gamma_ur`, which these methods call; they are listed under
[finding 012](012-statrs-gamma-cdf.md). The quantile is
[finding 036](036-statrs-chi-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Chi, ContinuousCDF};
use statrs::function::gamma::gamma_lr;

fn main() {
    let chi = |k| Chi::new(k).unwrap();
    println!("{:?}", gamma_lr(0.5, 1e-323)); // 0, should be 3.5470097772072545e-162
    println!("{:?}", gamma_lr(0.5, 2.5e-323)); // 0, should be 5.608314889430661e-162
    println!("{:?}", chi(1).cdf(4e-162)); // 0, should be 3.1915382432114613e-162
    println!("{:?}", chi(1).cdf(1e-160)); // should be 7.978845608028654e-161
    println!("{:?}", chi(1).cdf(8.696865935987788e-161)); // should be 6.939095057697016e-161
    println!("{:?}", std::panic::catch_unwind(|| chi(3).sf(1e155)).is_err()); // true: panics, should be 0
    println!("{:?}", std::panic::catch_unwind(|| chi(1).cdf(1e-170)).is_err()); // true: panics, should be 7.978845608028654e-171
    println!("{:?}", chi(10).cdf(1e-170)); // should be 0.0; panics
}
```

This is the program that was run (release and debug mode): it prints `0.0`,
`0.0`, `0.0`, `7.978801194389789e-161`, `6.941634974361731e-161`, `true`
and `true` (after the panic messages) and then panics with
``called `Result::unwrap()` on an `Err` value: XInvalid`` at
`src/function/gamma.rs:267`.

## Expected

| k | x | Method | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|---|
| 1 | 4e-162 | `cdf` | 0.0 | 3.191538243211461e-162 | 3.547009777207209e-162 |
| 1 | 1e-160 | `cdf` | 7.978801194389789e-161 | 7.978845608028654e-161 | 7.978801194389789e-161 |
| 1 | 8.696865935987788e-161 | `cdf` | 6.941634974361731e-161 | 6.939095057697016e-161 | (not on the grid) |
| 1 | 1e-170 | `cdf` | panic (XInvalid) | 7.978845608028654e-171 | 0.0 |
| every grid k | 5e-324, 1e-170 | `cdf` / `sf` | panic (XInvalid) | 0 / 1 (k = 1: 0.8 x / 1) | 0 / 1 |
| every grid k | 1e155, 1e300 | `cdf` / `sf` | panic (XInvalid) | 1 / 0 | 1 / 0 |

Second opinion: SciPy (`scipy.stats.chi`) gets the 52 points of group 1
right except k = 1, x = 1e-170 (0, as it forms x * x too), and has the same
problem as statrs in group 3 (the same value) and a similar one in group 2.
There the closed form for k = 1, F(x) = erf(x / sqrt 2) = x sqrt(2 / pi)
(1 - x^2 / 6 + ...), needs no square of x and gives the mpmath values; it is
tested against the dataset's route at x = 5e-324, 1e-170, 4e-162 and 1e-160
in `reference/tests/test_functions.py`.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route, P(k/2, x^2/2) and Q(k/2, x^2/2) with
   x^2/2 at 60 extra bits (`_mp_chi_cdf`, `_mp_chi_sf`, on the gamma_cdf
   and gamma_sf series and continued fraction), is tested against
   `mpmath.gammainc` at x^2/2, the closed forms for k = 1 (erf), k = 2
   (e^(-x^2/2)) and k = 3 (DLMF 8.4.6, computed from x itself, so also where
   x * x is subnormal or 0 in a double), the finite sum DLMF 8.4.10 at
   k = 1e4, and, at x = 0, 1e155 and 1e300, the end values, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 5,000 random points (Python `random.Random(20261016)`,
   shared with finding 036's scan): k = round(10^u) with u uniform in
   [0, 4]; x = 10^w with w uniform in [-175, -155] or [150, 160] (10%), else
   x = sqrt(2 t) with t = a + z sqrt(a), a = k/2, z uniform in [-40, 40]
   (50%; where that is not positive, t = a 10^v with v uniform in [-30, 0])
   or t = a 10^v with v uniform in [-30, 4] (40%); reference by the
   dataset's route at 30 digits. `cdf` and `sf` each panic at 333 points:
   170 with x <= 2.47e-162 and 163 with x >= 1.52e154, at k from 1 to 9889;
   no x between 2.79e-162 and 1.15e154 panics. `cdf` is also wrong at 3
   points with k = 1 and a normal result where x * x / 2 is subnormal
   (3.7e-4 at x = 8.7e-161, 1.4e-6 at 9.3e-160, 3.0e-8 at 6.4e-159). The
   other failures of the scan are finding 012's (section "Notes").
5. **Root cause** (from the source, `src/distribution/chi.rs` and
   `src/function/gamma.rs`): `cdf` returns
   `gamma::gamma_lr(k / 2.0, x * x / 2.0)` and `sf` returns
   `gamma::gamma_ur(k / 2.0, x * x / 2.0)` for every finite x > 0.
   - **Cause 1:** for x below 2.72e-162 the argument rounds to 0, and from
     1.34e154 up `x * x` is +inf. `checked_gamma_lr` and `checked_gamma_ur`
     reject x <= 0 and x = +inf with `Err(XInvalid)`, and the unchecked
     functions `unwrap` that error (finding 032's panic, reached from a
     positive x).
   - **Cause 2:** `checked_gamma_lr` returns 0 when
     `prec::ulps_eq!(x, 0.0, epsilon = 0.0)` holds, which is true for the
     five smallest positive doubles (up to 2.5e-323; measured: the first
     nonzero `gamma_lr(0.5, x)` is at x = 3e-323). P(1/2, x) there is about
     1.13 sqrt(x), 3.5e-162 at x = 1e-323, not 0. For k >= 2 the true
     value is below 2.5e-323 and 0 is right.
   - **Cause 3:** a subnormal `x * x / 2` carries fewer than 53 bits (about
     11 at x = 1e-160), and for k = 1 the result, about 0.8 x, is a normal
     number that inherits that error.

   Confirmed with a copy (scratch, not committed): returning 0 / 1 from
   `cdf` and 1 / 0 from `sf` where `x * x / 2` is 0 or +inf, and, for k = 1
   with x^2/2 below 1e-300, x sqrt(2 / pi) (whose next term is below 1e-300
   relative), leaves 237 of the scan's 10,000 evaluations failing, all of
   them finding 012's (k >= 1454 or subnormal results).

## Notes

`gamma_lr`'s treatment of the five smallest subnormals as 0 (cause 2) is
also there for direct callers and for `Gamma::cdf` with a subnormal x; no
dataset of this project has such an x for those. The off-grid scan's other
`cdf` / `sf` failures (112 and 128) are finding 012's: errors in normal
results for k from 1454 up (95 and 127 points, worst 7.3e-12 and 9.9e-12),
and subnormal results returned as 0 (17 and 1).
