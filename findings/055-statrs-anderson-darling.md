# 055: statrs anderson_darling takes the upper tail as 1 - cdf: wrong statistics, and inf with a NaN p-value, for values far in the upper tail

- **Crate / version:** statrs 0.19.1
  (`statrs::stats_tests::anderson_darling::anderson_darling`)
- **Severity:** high (silently wrong answers: a statistic 0.11% off and a
  p-value 5.6% off with one value at x = 36 in an exponential sample; inf
  and NaN from x = 37.43 on, where the true statistic is an ordinary number
  such as 8.98)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 48 of 128 points of `reference/data/anderson.jsonl`
  (tolerance 3e-12, for the statistic and the p-value), all exponential
  samples (dist 1) whose largest value is 20 or more. The exact list is
  `ANDERSON` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Every uniform
  sample passes (worst 1.4e-12), and so do the exponential samples whose
  largest value is 10 or less. One cause (section "Root cause"):
  **the upper tail is formed as `1 - cdf`**. The statistic needs
  log(1 - F(x)) at every value; statrs computes `1.0 - dist.cdf(x)`, which
  keeps only an absolute accuracy of about 1e-16. Where the true tail S(x)
  is small this loses digits: relative errors of the statistic from
  3.5e-12 (x = 20) through 4.8e-6 (x = 30) to 1.1e-3 (x = 36), and of the
  p-value up to 5.6e-2. Once F(x) rounds to 1 (x >= 37.43 for the
  exponential), `1 - cdf` is 0, its log is -inf, and statrs returns
  A^2 = inf with p-value NaN (30 points at x = 38 to 800).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Exp};
use statrs::stats_tests::anderson_darling::anderson_darling;

fn main() {
    let e = Exp::new(1.0).unwrap();
    for far in [20.0, 30.0, 36.0, 40.0] {
        let x = [far, 1.233642578125, 1.100341796875, 0.967041015625, 0.833740234375];
        println!("{far}: {:?}", anderson_darling(&x, &e).unwrap());
    }
    // The cause: 1 - cdf against sf.
    println!("{:e} {:e}", 1.0 - e.cdf(30.0), e.sf(30.0));
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below (statistic, p-value). The last
line prints `9.35918009759007e-14 9.357622968840175e-14`: `1 - cdf` is
already 1.7e-4 off at x = 30.

## Expected

The sample is 5 values of the grid's exponential pattern with the largest
replaced by `far`.

| far | Crate | mpmath (true value) | SciPy |
|---|---|---|---|
| 20 | 4.9826346520567135, 3.559552082310104e-15 | 4.982634648049917, 3.559552179211129e-15 | 4.982634648049915, 3.559552179211168e-15 |
| 30 | 6.982601366675059, 5.010049034900131e-21 | 6.982634644340009, 5.008935490292334e-21 | 6.982634644340008, 5.008935490292387e-21 |
| 36 | 8.191365322163273, 1.6269867393748335e-24 | 8.18263464433984, 1.7236550919921148e-24 | 8.182634644339839, 1.7236550919921537e-24 |
| 40 | inf, NaN | 8.982634644339841, 8.860533181152329e-27 | 8.98263464433984, 8.860533181152467e-27 |

Second opinion: SciPy's statistic for a fully specified distribution (the
one `scipy.stats.goodness_of_fit` uses, with `logcdf` and `logsf`) agrees
with mpmath within 2.1e-13 at every point of the dataset; SciPy has no
counterpart of the p-value, so the same published formula is evaluated in
doubles on SciPy's statistic (within 2.1e-13 too).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The statistic is the sum
   formula with log F and log S formed separately (log(-expm1(-x)) and -x
   for the exponential). New tests in `reference/tests/test_functions.py`
   check it against the defining integral n int (F_n - u)^2 / (u (1 - u)) du
   in closed form at every grid point (to 1e-45), against numerical
   quadrature of that integral, and against the one-value closed form; the
   p-value against the formula typed by hand from D'Agostino and Stephens
   (1986, Table 4.9), one point per branch.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random samples (NumPy `default_rng(20261103)`), n
   log-uniform in 1..300, judged at the larger of 3e-12 and each sample's
   own rounding bound (`anderson_rounding_bound`), true values from the
   reference at 30 and 50 digits (which agree at every point); 100 samples
   whose formula p-value exceeds the double range were left out. Half are
   uniform samples u^g (u uniform on (0, 1), g = exp(N(0, 0.6))): none
   fails. Half are exponential samples 1/2 + E with E exponential of a
   log-uniform scale in 0.1..20, and in 30% of them one value replaced by a
   log-uniform value in 5..1000: 570 fail, every one with a largest value of
   9.07 or more.
5. **Root cause** (from the source, `src/stats_tests/anderson_darling.rs`,
   confirmed with a copy that has the remedy):
   `f64::ln(1.0 - dist.cdf(f_obs[n - 1 - i]))`. Replacing it by
   `dist.sf(f_obs[n - 1 - i]).ln()` (statrs's `Exp::sf` is exact on the
   `expon_sf` grid) makes the copy pass every grid point and every random
   point with a largest value below 709, and gives the true values in the
   table above. Beyond about 709 even the remedy fails: e^-x is subnormal
   (fewer digits) or 0 in doubles, though log S(x) = -x is an ordinary
   number. The copy fails the 12 grid points at x = 745 and 800 (inf and
   NaN at 800) and 17 random points (x from 727 to 984). There, only a
   logarithm of the tail computed directly (SciPy's `logsf`) works, and
   statrs's `ContinuousCDF` trait has no such method; statrs fails those
   points already through `1 - cdf`, so they are listed under this cause.

## Notes

The p-value is D'Agostino and Stephens's (1986, Table 4.9) empirical
formula for a **normal** sample whose mean and variance were estimated
from it. statrs applies it to whatever distribution is passed, including
fully specified ones (its own tests use `Gamma::new(2.0, 1.0)`), for which
the published p-value is a different one (the distribution of A^2 itself,
Anderson and Darling 1954; Marsaglia and Marsaglia 2004). This dataset
checks only that the formula is computed as written, so this is not
counted as a failure. Neither are two flaws of the formula itself, which
statrs reproduces: its p-value is smallest at z = 153.47 (about 2e-190)
and **rises again** above that, so between z = 153.47 and 306.71 a worse
fit gets a larger p-value (z = 300 gives 5.7e-17; 9 grid points lie
there), and above z = 306.71 it exceeds 1 (one grid point: p = 1.0e60 at
z = 329).
A user reading the p-value of a test against a fully specified
distribution gets a number from the wrong reference distribution. The
statrs function has no documentation saying which case it covers.

Not on the grid: a NaN in the sample makes the sort's
`partial_cmp(...).unwrap()` panic (there is no `NaNPolicy` argument, unlike
`ttest_onesample`, `f_oneway`, `skewtest` and the KS tests); a value
outside the support (cdf 0 or 1) gives inf, the true A^2 then being
infinite too. Only the uniform and exponential
distributions are used; with others the statistic also carries that
distribution's own cdf error (for example finding 009 for `Normal`, or
finding 023 for `Exp` below x = 1/2).
