# 009: statrs Normal cdf and sf have about 10 correct digits

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Normal`, methods
  `cdf` and `sf` of `ContinuousCDF`; the free functions `cdf_unchecked` and
  `sf_unchecked` are the same code)
- **Severity:** low (inaccuracy, up to 1.03e-10 relative, about 10 correct digits
  instead of 15 to 16)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 155 of 441 points of `reference/data/norm_cdf.jsonl` and
  155 of 441 points of `reference/data/norm_sf.jsonl` (tolerance 4e-13); the
  exact lists are `NORMAL_CDF` and `NORMAL_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  `cdf` fails at every grid point with -37.5 <= x <= -0.75 (4.9e-11 to 1.0e-10)
  and 0.75 <= x <= 2.25 (the same absolute error, diluted by a value near 1:
  5.4e-13 at x = 2.25); `sf` is the mirror image. |x| <= 0.5 passes, and so do
  the subnormal and zero results beyond x = -37.5, where the error metric is
  relative to 2.2e-308. The quantile `inverse_cdf` is **not** affected: it
  passes all 138 points of `norm_ppf` (worst 2.6e-16).
  `statrs::stats_tests::skewtest::skewtest` takes its p-values from
  `Normal::cdf`, so it inherits this error: 100 of 210 points of
  `reference/data/skewtest.jsonl` (tolerance 2e-12) fail only through it
  (p-values 2.4e-11 to 1.03e-10 off; list `SKEWTEST_ERFC`; added
  2026-10-15). Its own errors are finding 029.
  `statrs::stats_tests::mannwhitneyu::mannwhitneyu` takes its asymptotic
  p-values from `Normal::cdf` too: 268 of 690 points of
  `reference/data/mannwhitneyu.jsonl` (tolerance 1e-13) fail only through it
  (p-values 1.1e-13 to 1.03e-10 off; list `MANNWHITNEYU_ERFC`; added
  2026-10-29). Its own errors are finding 050.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Normal};

fn main() {
    let n = Normal::new(0.0, 1.0).unwrap();
    println!("{:?}", n.cdf(-0.75)); // 0.22662735239943582, should be 0.2266273523768682
    println!("{:?}", n.sf(1.0)); // 0.15865525394505725, should be 0.15865525393145705
}
```

The program that was run (scratch, not committed) reads x values from
standard input and prints `cdf(x)` and `sf(x)` with `{:?}`; it gives the
"Crate" column below.

## Expected

| Input x | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| -0.75 | `cdf` | 0.22662735239943582 | 0.2266273523768682 | 0.2266273523768682 | 1.0e-10 |
| -10 | `cdf` | 7.619853024559641e-24 | 7.619853024160525e-24 | 7.61985302416047e-24 | 5.2e-11 |
| -37.5 | `cdf` | 4.605353009356527e-308 | 4.605353009581955e-308 | 4.605353009581954e-308 | 4.9e-11 |
| 2.25 | `cdf` | 0.9877755273454866 | 0.9877755273449553 | 0.9877755273449553 | 5.4e-13 |
| 1 | `sf` | 0.15865525394505725 | 0.15865525393145705 | 0.15865525393145707 | 8.6e-11 |
| 2.5 | `sf` | 0.006209665325512105 | 0.006209665325776135 | 0.006209665325776132 | 4.3e-11 |

SciPy (`scipy.stats.norm`) agrees with mpmath to 7.3e-15 or better at every
row (second opinion).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route Phi(x) = erfc(-x / sqrt 2) / 2 is tested
   against `mpmath.ncdf` and against Gamma(1/2, x^2/2) / (2 sqrt(pi)), and
   cdf + sf = 1, in `reference/tests/test_functions.py`.
3. **Second opinion**: SciPy within tolerance at every row of "Expected".
4. **Off-grid**: 20,000 x drawn uniformly from [-38.5, 38.5] (Python
   `random.Random(20261014).uniform`), each compared with mpmath at 50 digits.
   `cdf` fails 10,107 of them and `sf` 9,955 (tolerance 4e-13), worst 1.03e-10
   (cdf at x = -0.709, sf at x = 0.708); failures span x = -37.65 to 2.35 for
   `cdf` and -2.36 to 37.65 for `sf`, matching the grid.
5. **Root cause** (confirmed): `Normal::cdf` is `0.5 * erfc((mean - x) / (sd *
   sqrt 2))` and `sf` is `0.5 * erfc((x - mean) / (sd * sqrt 2))`
   (`src/distribution/normal.rs`, `cdf_unchecked` and `sf_unchecked`), so they
   inherit the error of statrs's `erfc`, [finding 004](004-statrs-erfc.md),
   for |x / sqrt 2| >= 0.5 (the intervals of `erf_impl` with the mistyped
   constants). Replacing the 13 interval constants of `erf_impl` by their
   float-rounded values (`0.3440242112_f32 as f64`, ..., as in Boost) in a copy
   of statrs 0.19.1 removes every failure: 0 of the 40,000 off-grid checks and
   0 of the 441 `norm_cdf` grid points fail, worst 2.4e-13 (the x / sqrt 2
   rounding that the tolerance allows for), and the same copy passes all 281
   `erfc` points (worst 3.8e-16). Fixing `erfc` fixes this finding.

## Notes

Only the standard normal (mean 0, sd 1) is checked; other parameters go
through the same `erfc` call after `(x - mean) / (sd * sqrt 2)`, which adds the
rounding of that quotient. A finding of its own (rather than a citation of
004) because users call `Normal::cdf`, not `erfc`, and the affected range and
size of the error differ.
