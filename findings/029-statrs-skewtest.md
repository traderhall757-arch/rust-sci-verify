# 029: statrs skewtest reports perfectly symmetric samples as skewed (z = 1.01) and returns small p-values as 0

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::skewtest::skewtest`)
- **Severity:** high (silently wrong answers: for the sample 1, 2, ..., 8, whose
  skewness is exactly 0, it returns z = 1.0108 and a two-sided p-value of
  0.312 instead of z = 0 and p = 1; a two-sided p-value of 6.6e-34 is
  returned as 0)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 49 of 210 points of `reference/data/skewtest.jsonl`
  (tolerance 2e-12, for the statistic and the p-value). The exact list is
  `SKEWTEST` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Three causes
  (section "Root cause"):
  1. **Zero skewness replaced by 1** (30 points: all 10 exactly symmetric
     samples, n = 8 to 100, each alternative): z is 1.011 to 1.040 instead of
     0, so the two-sided p-value is about 0.31 instead of 1 and the one-sided
     ones about 0.16 and 0.84 instead of 0.5. SciPy does the same (it is
     where statrs's code comes from).
  2. **The p-value as 1 - cdf** (16 points): two-sided and "greater"
     p-values below about 1e-5 lose digits, and below about 1e-16 are 0.
  3. **Plain summation** (3 points, one sample of n = 100 with z = 0.030):
     the statistic is 6.8e-12 off, as in finding 026.

  Another 100 points of the dataset fail only through statrs's normal
  cdf, whose `erfc` is off by up to 1.03e-10 (finding 009); they are listed
  under finding 009 as `SKEWTEST_ERFC`.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::skewtest::skewtest;
use statrs::stats_tests::{Alternative, NaNPolicy};

fn main() {
    let t = |x: Vec<f64>, a: Alternative| skewtest(x, a, NaNPolicy::Error).unwrap();
    let symmetric: Vec<f64> = (1..=8).map(f64::from).collect();
    println!("{:?}", t(symmetric.clone(), Alternative::TwoSided)); // should be (0.0, 1.0)
    println!("{:?}", t(symmetric, Alternative::Greater)); // should be (0.0, 0.5)
    let mut outlier = vec![0.0; 99];
    outlier.push(1.0);
    println!("{:?}", t(outlier.clone(), Alternative::TwoSided));
    println!("{:?}", t(outlier, Alternative::Less));
    let skewed: Vec<f64> = (1..=30).map(|i| f64::from(i * i)).collect();
    println!("{:?}", t(skewed.clone(), Alternative::Greater));
    println!("{:?}", t(skewed, Alternative::Less));
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| Sample, alternative | z, crate | z, mpmath | p, crate | p, mpmath (true value) | p, SciPy |
|---|---|---|---|---|---|
| 1..8, two-sided | 1.0108048609177787 | 0 | 0.3121098361688086 | 1 | 0.3121098361421897 |
| 1..8, greater | 1.0108048609177787 | 0 | 0.1560549180844043 | 0.5 | 0.15605491807109484 |
| 99 zeros and a 1, two-sided | 12.1382011360467 | 12.138201136046696 | 0.0 | 6.6266876328461014e-34 | 6.626687632845903e-34 |
| 99 zeros and a 1, less | 12.1382011360467 | 12.138201136046696 | 1.0 | 1.0 | 1.0 |
| 1, 4, 9, .., 900, greater | 1.5479267804011096 | 1.5479267804011108 | 0.06081996303311887 | 0.06081996303049563 | 0.06081996303049576 |
| 1, 4, 9, .., 900, less | 1.5479267804011096 | 1.5479267804011108 | 0.9391800369668811 | 0.9391800369695044 | 0.9391800369695043 |

The last two rows (4.3e-11 and 2.8e-12 off) are finding 009's `erfc`
error, shown for comparison. Second opinion: SciPy (`scipy.stats.skewtest`)
agrees with mpmath within 3.0e-14 at every row except the symmetric sample,
where it makes the same substitution. For that row the true values need no
computation: the sample skewness of 1..8 is exactly 0 (the cubed deviations
cancel in pairs), so z = delta asinh(0) = 0, and the normal distribution's
two-sided p-value at 0 is 1. The published example (SciPy's documentation,
also in statrs's tests: z = 2.7788579769903414, p = 0.005455036974740185)
is reproduced by the reference route (tested).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route (`_mp_skewtest`)
   follows D'Agostino, Belanger and D'Agostino (1990), eqs. 4-9, with z =
   delta asinh(y / alpha) (the same as SciPy's delta log(y / alpha +
   sqrt((y / alpha)^2 + 1))) and the p-value from erfc. New tests in
   `reference/tests/test_functions.py` check the skewness against exact
   rational arithmetic, the published example, the p-value forms against
   `mpmath.ncdf`, and that the symmetric samples have z = 0 with moments
   exact in double arithmetic.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random samples (Python `random.Random(20261018)`):
   n uniform in 8..100; 10% exactly symmetric (random integers 1..1000 and
   their negatives); the rest loc + scale (g + c g^2 / 4) with g standard
   normal, c uniform in [-2, 2], loc 0 or uniform in [-20, 20], scale 10^u
   with u uniform in [-2, 2], rounded to 6 significant digits; alternative
   uniform among the three; mpmath at 50 and 100 digits (they agreed at all
   points); 44 whose rounding bound (the tolerance reason's) exceeds the
   tolerance are left out. Of the 2,956 left: all 331 symmetric samples
   fail (z from 1.011 to 1.041, never 0); 293 more fail through cause 2
   (one returned as 0) and 11 through cause 3 (up to 6.1e-12); another 1,374
   only through finding 009's `erfc` (up to 1.2e-10).
5. **Root cause** (from the source, `src/stats_tests/skewtest.rs`, confirmed
   with a copy that is bit for bit equal to statrs on all 210 grid points):
   - **Cause 1:** `if y == 0.0 { y = 1.0; }` before
     `zscore = delta * (y / alpha + ((y / alpha).powi(2) + 1.0).sqrt()).ln()`.
     The transform is defined at y = 0 (it gives 0), so the substitution only
     changes the answer. SciPy has the same line (`y = np.where(y == 0, 1, y)`
     in `scipy.stats.skewtest`). Removing it makes every symmetric grid
     sample and all 331 random ones return z = 0, p = 1 (two-sided) and 0.5.
   - **Cause 2:** `2.0 * (1.0 - norm_dist.cdf(zscore.abs()))` and
     `1.0 - norm_dist.cdf(zscore)`, as in finding 026. Using
     `norm_dist.sf(..)` instead removes this cause; the p-values that then
     still fail (all within 1.04e-10) are exactly the `erfc` error of
     finding 009, which `Normal::sf` shares.
   - **Cause 3:** the mean in `calc_root_b1` is a plain sum; when the
     skewness is near 0 the cubed deviations cancel and that rounding
     becomes relative error in z.

## Notes

Cause 1 matters most in practice: textbook and test data are often
perfectly symmetric (equally spaced values, mirrored values), and the test
then reports skewness (one-sided p = 0.16) where there is none. The same
code is used by every caller of `skewtest`.
