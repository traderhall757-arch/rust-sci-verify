# 026: statrs ttest_onesample returns p-values below about 1e-16 as 0 and loses digits below about 4e-5

- **Crate / version:** statrs 0.19.1
  (`statrs::stats_tests::ttest_onesample::ttest_onesample`)
- **Severity:** high (a silently wrong answer: a two-sided p-value of
  1.8e-24 is returned as exactly 0, and every p-value below about 4e-5 for
  the two-sided and "greater" tests loses digits: 2.9% off at 1.8e-15,
  3.5e-7 at 7.9e-11)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 60 of 324 points of `reference/data/ttest_1samp.jsonl`
  (tolerance 5e-12, for the statistic and the p-value): 42 two-sided and 18
  "greater" tests: every grid point of those two alternatives whose p-value
  is below 8.7e-6 (the next one up, 2.2e-5, passes). The exact list is `TTEST_1SAMP` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  - **p below about 1e-16** (34 grid points): the p-value is exactly 0.
  - **p from about 1e-16 to about 4e-5** (26 grid points): the p-value has
    an absolute error of about 1.1e-16 (two-sided 2.2e-16), so a relative
    error of about 1e-16 / p.
  - The "less" alternative is correct at every grid point, and the statistic
    t is within tolerance at every grid point.
  - Off the grid only, a second, small cause: when the sample mean nearly
    equals mu, the statistic loses up to 2.3e-11 (plain summation of the
    sample; section "Root cause").

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::ttest_onesample::ttest_onesample;
use statrs::stats_tests::{Alternative, NaNPolicy};

fn main() {
    let x = vec![-0.5, 0.419, 0.338, 0.257, 0.176, 0.095, 0.014, -0.067, -0.148, -0.229];
    let t = |mu: f64, alt: Alternative| ttest_onesample(x.clone(), mu, alt, NaNPolicy::Error).unwrap();
    println!("{:?}", t(100.0, Alternative::TwoSided)); // p should be 1.782747872274205e-24
    println!("{:?}", t(100.0, Alternative::Less)); // p = 8.913739361371025e-25 (right)
    println!("{:?}", t(-30.0, Alternative::Greater)); // p should be 4.465095086022772e-20
    println!("{:?}", t(10.0, Alternative::TwoSided)); // p should be 1.829536608912975e-15
    println!("{:?}", t(-3.0, Alternative::TwoSided)); // p should be 7.869746574097906e-11
    println!("{:?}", t(1.0, Alternative::TwoSided)); // p should be 1.8190622058847864e-6
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below. The sample is the grid's n = 10,
loc = 0 sample.

## Expected

| mu | Alternative | t (crate = mpmath) | p, crate | p, mpmath (true value) | p, SciPy | Rel. error of p |
|---|---|---|---|---|---|---|
| 100 | two-sided | -1123.675509800379 | 0.0 | 1.782747872274205e-24 | 1.782747872274206e-24 | 1 |
| 100 | less | -1123.675509800379 | 8.913739361371047e-25 | 8.913739361371025e-25 | 8.91373936137103e-25 | 2.5e-15 |
| -30 | greater | 337.62141334783126 | 0.0 | 4.465095086022772e-20 | 4.4650950860227785e-20 | 1 |
| 10 | two-sided | -112.0084091593103 | 1.7763568394002505e-15 | 1.829536608912975e-15 | 1.8295366089129753e-15 | 2.9e-2 |
| -3 | two-sided | 34.1212831555107 | 7.869749296673945e-11 | 7.869746574097906e-11 | 7.869746574097918e-11 | 3.5e-7 |
| 1 | two-sided | -10.841699095203452 | 1.819062205976607e-6 | 1.8190622058847864e-6 | 1.819062205884787e-6 | 5.0e-11 |

(The crate's t at mu = -3 and 1 is one ulp from mpmath's, within
tolerance.) Second opinion: SciPy (`scipy.stats.ttest_1samp`) agrees with
mpmath within 1.6e-15 at every row, and the "less" p-value of the same t
from statrs itself is right. For the statistic, t^2 is also computed in
exact rational arithmetic from the doubles (tested).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route
   (`_mp_ttest_1samp`) forms the mean, variance and t with 60 extra bits and
   takes the p-value from the Student's t routes of the `t_cdf` / `t_sf`
   datasets (already tested against `mpmath.betainc`, the series and closed
   forms). New tests in `reference/tests/test_functions.py` check t^2 against
   exact rational arithmetic (`fractions.Fraction`) and the p-value against
   the closed forms for 1 and 2 degrees of freedom (n = 2 and 3) at every
   such grid point, all three alternatives.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random samples (Python `random.Random(20261016)`):
   n uniform in 2..200; x_i = loc + scale g_i with g_i standard normal,
   loc uniform in [-10, 10], scale 10^u with u uniform in [-3, 3]; mu chosen
   so that t is about +-10^v with v uniform in [-2, 3.5]; alternative
   uniform among the three; all rounded to 6 significant digits; mpmath at
   50 and 100 digits (they agreed at all points). 411 samples whose
   rounding bound (the tolerance reason's) exceeds the tolerance are left
   out, which leaves 2,589. Two-sided and "greater" p-values fail at 618
   points, 488 of them returned as 0 (every true p below 3.8e-17 is 0);
   every p from 1e-16 to 1e-6 fails (112 of 112) and 17 of 37 from 1e-6 to
   1e-4; worst errors 0.91, 2.0e-5, 4.3e-9 and 6.9e-11 per range
   [1e-16, 1e-12), [1e-12, 1e-8), [1e-8, 1e-6), [1e-6, 1e-4). No "less"
   p-value fails. The statistic fails at 16 points (up to 2.3e-11), all with
   the sample mean within about 1e-4 relative of mu (kappa_t from 2.2e4 to
   4.4e4).
5. **Root cause** (from the source, `src/stats_tests/ttest_onesample.rs`,
   confirmed by experiment):
   - **Cause 1:** the p-value is `2.0 * (1.0 - t_dist.cdf(tstat.abs()))`
     for the two-sided test and `1.0 - t_dist.cdf(tstat)` for "greater".
     The cdf near 1 is rounded to a double with an absolute error of up to
     1.1e-16, which the subtraction keeps; below about 1.1e-16 the cdf
     rounds to exactly 1 and the p-value is 0. The "less" test uses
     `t_dist.cdf(tstat)` directly, which is accurate in the lower tail.
     Using `t_dist.sf(..)` instead of `1.0 - t_dist.cdf(..)` (statrs's
     `StudentsT::sf`) removes all 60 grid failures and 616 of the 618
     off-grid p-value failures.
   - **Cause 2 (off the grid):** the mean is a plain left-to-right sum, with
     an absolute error of up to about n 2^-53 times the mean of |x_i|; when
     the mean nearly equals mu, mean - mu cancels and that error becomes
     relative in t (n kappa_t 2^-53 instead of a careful method's
     (kappa_t + n) 2^-53). It matters only for samples of over about 40
     points whose mean is within about 1e-4 of mu, where t is small and the
     p-value near 1/2 or 1, so it is of little practical consequence; it
     moves the p-value beyond tolerance at 2 scan points (1.0e-11 at
     p = 0.38, and 9.8e-12 at p = 1.9e-9 once cause 1 is fixed). A compensated
     (Neumaier) sum for the mean, with the rest unchanged, together with
     cause 1's fix, leaves 0 of the 2,589 scan points failing (statistic and
     p-value).

## Notes

The other statrs hypothesis tests compute their p-values the same way:
`chisquare` as `1.0 - chi_dist.cdf(stat)` and `f_oneway` as
`1.0 - f_dist.cdf(fstat)` (to be checked). SciPy computes the one-sample
p-value from the survival function (`scipy.special.stdtr` of -|t|).
