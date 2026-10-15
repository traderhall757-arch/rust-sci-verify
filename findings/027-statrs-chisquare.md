# 027: statrs chisquare returns p-values below about 1e-16 as 0 and loses digits below about 1e-5

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::chisquare::chisquare`)
- **Severity:** high (a silently wrong answer: a p-value of 2.8e-43 is
  returned as exactly 0, and every p-value below about 1e-5 loses digits:
  `chisquare(&[1000, 1300], None, None)` is 1.0e-8 off, `&[100, 200]` 1.7e-9)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 29 of 175 points of `reference/data/chisquare.jsonl`
  (tolerance 5e-12, for the statistic and the p-value) with `f_exp` given, and
  14 of the 76 points with equal expected frequencies when run with
  `f_exp = None` (the default). The exact lists are `CHISQUARE` and
  `CHISQUARE_EQUAL` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  - **p below about 1e-16** (22 of the 29): the p-value is exactly 0.
  - **p from about 1e-16 to about 1e-5** (7 of the 29): an absolute error of
    up to about 1.1e-16, so a relative error of about 1e-16 / p.
  - The statistic is within tolerance at every grid point and every random
    point.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::chisquare::chisquare;

fn main() {
    println!("{:?}", chisquare(&[492500, 506285], None, None)); // p should be 2.7935024579296593e-43
    println!("{:?}", chisquare(&[1000, 1300], None, None)); // p should be 3.9641597078900157e-10
    println!("{:?}", chisquare(&[100, 200], None, None)); // p should be 7.764036537930677e-9
    println!("{:?}", chisquare(&[16, 18, 16, 14, 12, 12], None, None)); // p = 0.8491450360846097 (right)
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below. The first sample is a grid
point; the last is the example in statrs's own tests.

## Expected

| f_obs (f_exp = None) | Statistic (crate = mpmath) | p, crate | p, mpmath (true value) | p, SciPy | Rel. error of p |
|---|---|---|---|---|---|
| 492500, 506285 | 190.2573877260872 | 0.0 | 2.7935024579296593e-43 | 2.793502457929687e-43 | 1 |
| 1000, 1300 | 39.130434782608695 | 3.9641601112805347e-10 | 3.9641597078900157e-10 | 3.9641597078900225e-10 | 1.0e-8 |
| 100, 200 | 33.333333333333336 | 7.764036524804396e-9 | 7.764036537930677e-9 | 7.764036537930667e-9 | 1.7e-9 |
| 16, 18, 16, 14, 12, 12 | 2.0 | 0.8491450360846101 | 0.8491450360846097 | 0.8491450360846096 | 4.6e-16 |

Second opinion: SciPy (`scipy.stats.chisquare`) agrees with mpmath within
1.0e-14 at every row, and at every point of the dataset within 1.9e-13.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route (`_mp_chisquare`)
   sums (o - e)^2 / e with 60 extra bits and takes the p-value
   Q(dof / 2, statistic / 2) from the gamma_sf route (the series and
   continued fraction tested against `mpmath.gammainc`, finite sums and
   closed forms). New tests in `reference/tests/test_functions.py` check the
   statistic against exact rational arithmetic and the p-value against the
   closed forms for 1 and 2 degrees of freedom (erfc(sqrt(x / 2)) and
   e^(-x / 2)) and `mpmath.gammainc`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random tests (Python `random.Random(20261017)`):
   k uniform in 2..50 categories with weights uniform in [0.5, 2] (3
   digits), total 10^u with u uniform in [1, 6], observed counts
   round(e (1 + delta g)) with g standard normal and delta 10^v, v uniform
   in [-3, -0.5]; expected frequencies the weights scaled to the counts'
   sum; ddof 0 or 1. 21 had all counts at their expectation or zero and were
   skipped; mpmath at 50 and 100 digits agreed at all 2,979 others; 5 whose
   rounding bound (the tolerance reason's) exceeds the tolerance are left
   out. Of the 2,974 left, 435 p-values fail (none of the statistics): 305
   are returned as 0 (every true p below 4.9e-17 is), and every failure has
   p at most 7.5e-6.
5. **Root cause** (from the source, `src/stats_tests/chisquare.rs`,
   confirmed by experiment): the p-value is `1.0 - chi_dist.cdf(stat)`. Near
   1 the cdf is rounded to a double with an absolute error of up to 1.1e-16,
   which the subtraction keeps; below about 1.1e-16 the cdf rounds to 1 and
   the p-value is 0. Using `chi_dist.sf(stat)` (statrs's `ChiSquared::sf`)
   with the same statistic passes every grid point and all but 2 of the
   2,974 random points; those 2 have subnormal true p-values (6.4e-317 and
   3.6e-315) that `sf` itself returns as 0, finding 012's flush.

## Notes

The same pattern is in `ttest_onesample` (finding 026) and `f_oneway`
(finding 028). SciPy computes this p-value with the survival function
(`scipy.special.chdtrc`).
