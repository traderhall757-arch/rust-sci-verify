# 013: statrs Gamma and ChiSquared inverse_cdf return about 1e-9 for small quantiles and lose digits near p = 1

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Gamma` and
  `statrs::distribution::ChiSquared`, method `inverse_cdf` of
  `ContinuousCDF`; `ChiSquared` and `Erlang` call `Gamma`'s)
- **Severity:** high (silently wrong answers: about 9.3e-10 instead of
  quantiles down to 1e-300, wrong by factors up to 4e298; 0.2% off near
  p = 1; 2% to 8% off in the far lower tail for large shape)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 170 of 432 points of `reference/data/gamma_ppf.jsonl`
  for `Gamma` (tolerance 2e-13), and 170 for `ChiSquared` with df = 2a (half
  its quantile is compared); the exact lists are `GAMMA_PPF_WHOLE` (111
  points with a whole-number shape) plus `GAMMA_PPF_FRACTIONAL` (59) and
  `CHI2_PPF` in `crates/check-statrs/src/expected.rs` (the `Gamma` and
  `ChiSquared` lists differ at 28 points, see cause 1), the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Only rate 1 is
  checked. `Erlang::inverse_cdf` calls `Gamma`'s, so it fails at exactly the
  111 points of `GAMMA_PPF_WHOLE` (of the 324 whole-number-shape points;
  added 2026-10-15). Five groups for `Gamma`, by symptom (`ChiSquared` in brackets):
  1. **Tiny quantiles** (65 [65] points; causes 1 and 4): wherever the true
     quantile is below about 1e-9 (small shape, or small p) the result is
     about 9.3e-10 (2^-30), e.g. 9.3e-10 for the quantile 1e-50 of a = 1.
     The error reaches 4e298.
  2. **p near 1** (50 [49] points; cause 3): 2e-8 at p = 1 - 1e-10 and 0.2%
     at p = 1 - 1e-15 for every shape; smaller errors from p = 0.9 on.
  3. **Far lower tail, shape >= 20** (15 [15] points, p <= 1e-50; causes 2
     and 4): up to 8% off for a >= 1000 (6934 instead of 6738 for a = 1e4,
     p = 1e-300), wrong by factors up to 2e4 for a = 20 to 100.
  4. **Bulk, small shape** (12 [15] points; cause 1): errors up to 4.8e-9
     (a = 0.5, p = 0.4) [2.3e-8].
  5. **Bulk, large shape** (28 [26] points, a >= 1e4; cause 5): up to
     1.3e-12.

`cdf` and `sf` have a smaller, separate problem:
[finding 012](012-statrs-gamma-cdf.md), which also causes groups 3 and 5.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ChiSquared, ContinuousCDF, Gamma};

fn main() {
    let g = |a| Gamma::new(a, 1.0).unwrap();
    println!("{:?}", g(1.0).inverse_cdf(1e-50)); // should be 1e-50
    println!("{:?}", g(0.5).inverse_cdf(1e-10)); // should be 7.853981633974484e-21
    println!("{:?}", g(0.01).inverse_cdf(0.5)); // should be 4.465535018910355e-31
    println!("{:?}", g(1.0).inverse_cdf(0.5000000001)); // should be 0.6931471807599453
    println!("{:?}", g(0.5).inverse_cdf(0.4)); // should be 0.13749794886422803
    println!("{:?}", g(1.0).inverse_cdf(0.999999999999999)); // should be 34.53957599234088
    println!("{:?}", g(3.0).inverse_cdf(0.9999999999)); // should be 29.145901390216682
    println!("{:?}", g(1e4).inverse_cdf(1e-300)); // should be 6737.687191590329
    println!("{:?}", g(1e3).inverse_cdf(5e-324)); // should be 218.26424714135766
    println!("{:?}", g(1e5).inverse_cdf(0.5)); // should be 99999.6666668642
    println!("{:?}", ChiSquared::new(1.0).unwrap().inverse_cdf(1e-10)); // should be 1.5707963267948968e-20
    println!("{:?}", ChiSquared::new(2.0).unwrap().inverse_cdf(0.999999999999999)); // should be 69.07915198468176
}
```

This is the program that was run (release mode); it gives the "Crate"
column below.

## Expected

| a | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| 1 | 1e-50 | 9.313225746154785e-10 | 1e-50 | 1.0000000000000025e-50 | 9.3e40 |
| 0.5 | 1e-10 | 9.313225746222546e-10 | 7.853981633974484e-21 | 7.85398163397448e-21 | 1.2e11 |
| 0.01 | 0.5 | 9.313225746154785e-10 | 4.465535018910355e-31 | 4.465535018910562e-31 | 2.1e21 |
| 1 | 0.5000000001 | 0.6931471801634228 | 0.6931471807599453 | 0.6931471807599453 | 8.6e-10 |
| 0.5 | 0.4 | 0.13749794820591593 | 0.13749794886422803 | 0.13749794886422795 | 4.8e-9 |
| 1 | 1 - 1e-15 | 34.485508770746364 | 34.53957599234088 | 34.53957599234088 | 1.6e-3 |
| 3 | 1 - 1e-10 | 29.145900796061394 | 29.145901390216682 | 29.145901390216682 | 2.0e-8 |
| 1e4 | 1e-300 | 6933.933702274718 | 6737.687191590329 | 6737.687191590329 | 2.9e-2 |
| 1000 | 5e-324 | 224.12500000031253 | 218.26424714135766 | 218.2643501534942 | 2.7e-2 |
| 1e5 | 0.5 | 99999.66666698278 | 99999.6666668642 | 99999.6666668642 | 1.2e-12 |
| `ChiSquared` df = 1 | 1e-10 | 9.313225746290308e-10 | 1.5707963267948968e-20 | | 5.9e10 |
| `ChiSquared` df = 2 | 1 - 1e-15 | 68.97101754230269 | 69.07915198468176 | | 1.6e-3 |

(1 - 1e-15 and 1 - 1e-10 are the doubles nearest those numbers. The
`ChiSquared` rows give the chi-squared quantile itself, twice the gamma
quantile at a = df / 2.)

Second opinion: SciPy (`scipy.stats.gamma`) agrees with mpmath within
4.7e-14 at every row except a = 1000, p = 5e-324, where it is also wrong
(4.7e-7; it fails all six p = 5e-324 points with a >= 100). For a = 1 and
a = 1/2 the quantile has closed forms, x = -log(1 - p) and x = erfinv(p)^2
(tested at p = 1e-300 to 1 - 1e-15), and for a = 1000, p = 5e-324 the
dataset's route is tested by P(a, x) = p to 40 digits
(`reference/tests/test_functions.py`).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The quantile is found by bisection and
   Newton's method on log P(a, e^u) = log p (log Q(a, e^u) = log(1 - p) for
   p > 1/2), on the cdf route of finding 012; tests check P(a, ppf(p)) = p
   (or Q = 1 - p) at a = 0.01 to 1e6 and p = 5e-324 to 1 - 1e-15, and the
   closed forms for a = 1 and 1/2.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points, drawn after the 20,000 of finding 012
   with the same generator (`random.Random(20261014)`): a = 10^u with u
   uniform in [-2, 6]; p uniform in (0, 1) (40%), p = 10^v with v uniform
   in [-300, -1] (30%), or p = 1 - 10^v with v uniform in [-15.5, -1]
   (30%); reference by the dataset's route at 30 digits.
   `Gamma::inverse_cdf` fails 1,613 (54%) and `ChiSquared::inverse_cdf`
   1,623; the worst is 9.3e-10 [4.7e-10] where the true quantile is below
   the smallest double (error 4.2e298). 431 have a true quantile below
   1e-9 (group 1); every call returned. With the fixes for causes 1 to 3
   below, 1,457 [1,451] pass; the 156 [172] left all have a >= 1.1e4
   [1.7e4] and errors up to 2.0e-12 [2.4e-12] (cause 5; statrs's cdf error
   at the true quantile, divided by x f(x), predicts the observed error to
   a median factor of 1.2 [1.3]).
5. **Root cause** (confirmed by experiment): `Gamma::inverse_cdf`
   (`src/distribution/gamma.rs`) brackets the quantile by halving 1 and
   doubling 2, then runs at most 100 steps of a safeguarded Newton method on
   `cdf(x) - p`, stopping when `prec::convergence` says two iterates are
   equal. A copy of the function (scratch, not committed; it reproduces
   statrs bit for bit on all 800 grid evaluations of both distributions)
   with switchable fixes shows:
   - **Cause 1, absolute stopping test** (groups 1 and 4):
     `prec::convergence` is `relative_eq!` with statrs's defaults
     `epsilon = 1e-9` (absolute) and `max_relative = 1e-14`, and `approx`
     accepts when either holds. So the iteration stops as soon as two
     iterates differ by less than 1e-9: quantiles below that come back as
     the bracket point reached by then (about 2^-30), and elsewhere the
     last step, up to 1e-9 in absolute size, is the error. `ChiSquared` runs
     the same search on a scale twice as large, so the absolute test cuts at
     different points: the two lists differ at 28 points. Stopping on a
     relative change of 4 ulps instead fixes 37 [38] of the 170 points by
     itself, and 88 [88] together with cause 2's fix.
   - **Cause 2, iteration cap** (group 3): in the far lower tail for large a
     the Newton step (cdf - p) / pdf is about x / (a - x), a few units, while
     the root is hundreds of units away, and the loop ends after 100 steps
     and returns the unconverged iterate without a warning. Allowing
     100,000 steps moves, e.g., a = 1e4, p = 1e-300 from 6934 to 6737.69 (all
     digits right).
   - **Cause 3, cdf - p near p = 1** (group 2): the search solves
     cdf(x) = p, but near 1 the cdf carries an absolute rounding error of
     about 1e-16, which moves the root by 1e-16 / f(x). In the upper tail
     f(x) is 1 - p times a modest factor (exactly 1 - p for a = 1), so the
     relative error is about 1e-16 / ((1 - p) x): 3e-3 for a = 1 at
     1 - p = 1e-15 (observed 1.6e-3). Solving sf(x) = 1 - p instead (1 - p
     is exact) fixes 49 [48] points by itself; with causes 1 and 2 fixed as
     well, 130 [130] of the 170 pass.
   - **Cause 4, subnormal p** (13 of the 40 left, all p = 5e-324): `cdf`
     returns 0 below about 1e-308 (finding 012, cause 2), so the search
     cannot see p. Without that cut-off it gets closer, but a double cdf
     near 5e-324 has a single significant bit, so only a search on log P
     (as the reference does) finds these quantiles.
   - **Cause 5, cdf accuracy** (group 5, the other 27 of the 40): the search
     returns a root of statrs's own `cdf`, whose relative error for
     a >= 1e4 is up to 1.2e-9 (finding 012, cause 1); divided by the
     condition number x f(x) / P(x) (hundreds near the mean) it gives the
     2e-13 to 1.3e-12 seen. Evaluating statrs's cdf error at the true
     quantile and dividing by x f(x) predicts errors of the same size
     (3.6e-14 to 3.4e-12) at these 27 points.
   So every failure has one of these causes: after the fixes for causes 1
   to 3, the 40 [40] points left are the 13 [13] at p = 5e-324 (cause 4)
   and 27 [27] with a >= 1e4 and errors of 2.2e-13 to 1.3e-12 (cause 5).

## Notes

statrs's `StudentsT::inverse_cdf` (finding 011) has an independent
implementation with similar weaknesses (absolute stopping, no iteration
limit). Boost.Math's `gamma_p_inv` starts from DiDonato and Morris's
approximation and uses Halley steps with a relative tolerance, choosing
P or Q by p; SciPy's `gammaincinv` likewise.
