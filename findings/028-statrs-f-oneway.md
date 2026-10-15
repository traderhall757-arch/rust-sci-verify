# 028: statrs f_oneway loses digits in F when the data are large compared with their spread, and returns small p-values as 0

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::f_oneway::f_oneway`)
- **Severity:** high (silently wrong answers: for two groups of three values
  near 1e5 with spread 0.1 the statistic is 0.80014 instead of 0.80000, and
  a p-value of 8.6e-15 is 0.28% off; p-values below about 1e-16 are 0)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 34 of 48 points of `reference/data/f_oneway.jsonl`
  (tolerance 5e-12, for the statistic and the p-value). The exact list is
  `F_ONEWAY` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes
  (section "Root cause"):
  1. **The statistic** (25 points, up to 8.5e-9): all 24 grid points with
     values around 100 (spread about 0.3 within groups), and 1 around 0
     where the groups are 30 apart (F = 1.2e5, 2.2e-11). The p-value
     inherits the error, even near 1 (p = 0.72).
  2. **The p-value as 1 - cdf** (9 more points with a correct statistic): 0
     for every p below about 1e-16 (10 grid points in all), digits lost
     below about 1e-5.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::f_oneway::f_oneway;
use statrs::stats_tests::NaNPolicy;

fn main() {
    let t = |g: Vec<Vec<f64>>| f_oneway(g, NaNPolicy::Error).unwrap();
    println!("{:?}", t(vec![vec![1000.1, 1000.2, 1000.4], vec![1000.2, 1000.3, 1000.6]]));
    println!("{:?}", t(vec![vec![100000.1, 100000.2, 100000.4], vec![100000.2, 100000.3, 100000.6]]));
    println!("{:?}", t(vec![vec![1.0, 2.0, 3.0], vec![101.0, 102.0, 103.0]]));
    println!("{:?}", t(vec![vec![1.0, 2.0, 3.0, 4.0], vec![101.0, 102.0, 103.0, 104.0], vec![51.0, 52.0, 53.0, 54.0]]));
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| Groups | F, crate | F, mpmath (true value) | p, crate | p, mpmath | Rel. error (F; p) |
|---|---|---|---|---|---|
| {1000.1, 1000.2, 1000.4}, {1000.2, 1000.3, 1000.6} | 0.7999999329447749 | 0.8 | 0.4216482729985762 | 0.4216482551761941 | 8.4e-8; 4.2e-8 |
| the same + 99000 | 0.8001373233392459 | 0.8000000000349246 | 0.4216117593375164 | 0.42164825516691157 | 1.7e-4; 8.7e-5 |
| {1, 2, 3}, {101, 102, 103} | 15000.0 | 15000.0 | 2.6654818907623223e-8 | 2.665481896163596e-8 | 0; 2.0e-9 |
| {1, 2, 3, 4}, {101, .., 104}, {51, .., 54} | 6000.0 | 6000.0 | 8.659739592076221e-15 | 8.635957388336973e-15 | 0; 2.8e-3 |

(The decimal inputs are rounded to doubles, so the true F of the first two
rows is not exactly 0.8.) Second opinion: SciPy (`scipy.stats.f_oneway`)
agrees with mpmath within 9.5e-16 at every row. On the grid SciPy is off at
4 points (it subtracts sums of squares too, after removing the overall
mean): 5.1e-12 in one statistic with F = 1.2e5 and 2.0e-11 to 7.5e-11 in
p-values below 1e-52; there the p-value is confirmed by `mpmath.betainc`
(tested).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route (`_mp_f_oneway`)
   forms the group means, the between- and within-group sums of squared
   deviations and F with 60 extra bits, and the p-value from the f_sf route
   (tested against `mpmath.betainc`, the series and closed forms). New tests
   in `reference/tests/test_functions.py` check F against exact rational
   arithmetic and the p-value against `mpmath.betainc` and, for two groups,
   against the two-sample t-test's p-value (F = t^2).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 2,000 random analyses (Python `random.Random(20261017)`,
   after finding 027's 3,000 tests): k uniform in 2..6 groups of 2..40
   values; values loc + shift_j + sd g (g standard normal), with loc 0 or
   +-10^u (u uniform in [-1, 3]), sd 10^v (v uniform in [-1, 1]), group
   shifts sd 10^w g (w uniform in [-3, 2]); all rounded to 6 significant
   digits; mpmath at 50 and 100 digits (they agreed at all points); 41 whose
   rounding bound (the tolerance reason's) exceeds the tolerance are left
   out. Of the 1,959 left:
   - the statistic fails at 322 (worst 4.5e-7). By |mean| / sd (sd within
     groups): 3 of 1,320 below 10, 105 of 413 from 10 to 100, 172 of 184
     from 100 to 1000, all 42 above 1000 (worst per range 5.4e-10,
     4.1e-8, 1.9e-7, 4.5e-7);
   - the p-value fails at 1,006, 623 of them returned as 0 (every true p
     below 4.9e-17 is). 696 of the 1,006 have a correct statistic.
5. **Root cause** (from the source, `src/stats_tests/f_oneway.rs`, confirmed
   by experiment):
   - **Cause 1:** the sums of squares are the textbook one-pass formulas
     `sst = tsq - g^2 / n` and `sse = ysq - tsq`, with `ysq` the sum of all
     x^2 and `tsq` the sum of (group sum)^2 / n_j. Both subtract large,
     nearly equal numbers when the values are large compared with their
     spread within groups (through an offset or very different group
     means), losing about (mean / sd)^2 2^-53 or more.
   - **Cause 2:** the p-value is `1.0 - f_dist.cdf(fstat)`, as in
     findings 026 and 027.
   - A copy with two-pass sums of squares (group means first, then squared
     deviations) and `FisherSnedecor::sf` passes all but 2 grid points and
     all but 48 of the 1,959 random points. 43 of those 48 and both grid
     points are p-values below 1e-59 off by up to 2.2e-11 (grid: 9e-12) with
     a correct statistic: statrs's F survival function itself, finding 021's
     prefactor error. The other 5 come from the copy's plain sums: 3
     statistics 1.02 to 1.18 times the tolerance (5.1e-12 to 5.9e-12) and 2
     p-values (6.0e-12 and 6.9e-12) from statistics just inside it.

## Notes

Even moderately offset data are affected: values around 1000 with a spread
of 0.1 (the first row) lose half the digits of F. SciPy subtracts the
overall mean before forming the sums of squares, which removes the offset
part of cause 1.
