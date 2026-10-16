# 038: statrs Weibull inverse_cdf gives 0, inf or NaN where scale^-shape leaves the double range

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Weibull`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with shape 3.5 and scale
  2^500, the median is inf instead of 2.95e150, and with shape 100 and
  scale 1/1024, `inverse_cdf(1e-50)` is 0 instead of 3.09e-4)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 214 of 1215 points of `reference/data/weibull_ppf.jsonl`
  (tolerance 1e-13; 45 of the 1215 are p = 1, tagged pole). The exact list
  is `WEIBULL_PPF` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Every failure
  has an intermediate outside the normal double range (2.2e-308 to
  1.8e308), although the answer is not:
  1. **scale^-shape overflows to inf** (100 points; scale 2^-500 with
     shape >= 2.1): 0 at every p.
  2. **scale^-shape underflows to 0** (78 points; scale 2^500 with
     shape >= 3.5): inf at every p > 0, NaN at p = 0 (0 / 0).
  3. **scale^-shape is subnormal** (19 points; scale 2^500, shape 2.1:
     2^-1050): inf, because L = -ln(1 - p) divided by it overflows.
  4. **L / scale^-shape underflows to 0 or a subnormal** (17 points; shape
     100 with scale 2^-10 and shape 2 with scale 2^-500, small p): 0, or up
     to 2.8e-4 off.

  Ordinary-looking parameters are affected too: with scale 2^-10 and shape
  100, every p below 3e-13 gives a wrong answer (0 below 2.5e-23; a scan
  of 3,101 p from 1e-320 to 1e-10 also found a few failures up to 1e-12). In the random scan
  below with scales from 2^-20 to 2^20, 16% of points fail.

`Weibull::cdf` and `sf` fail for the same reason:
[finding 037](037-statrs-weibull-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Weibull};

fn main() {
    let tiny = 2f64.powi(-500); // 3.054936363499605e-151
    let huge = 2f64.powi(500); // 3.273390607896142e150
    let w = Weibull::new(3.5, tiny).unwrap(); // scale^-3.5 = 2^1750 overflows
    println!("{:?}", w.inverse_cdf(0.5)); // should be 2.7512099172225064e-151
    let w = Weibull::new(3.5, huge).unwrap(); // scale^-3.5 = 2^-1750 underflows to 0
    println!("{:?}", w.inverse_cdf(0.5)); // should be 2.9479451064801336e150
    println!("{:?}", w.inverse_cdf(0.0)); // should be 0.0
    let w = Weibull::new(2.1, huge).unwrap(); // scale^-2.1 = 2^-1050 is subnormal
    println!("{:?}", w.inverse_cdf(1e-5)); // should be 1.3615904136475406e148
    let w = Weibull::new(100.0, 2f64.powi(-10)).unwrap(); // -ln(1 - p) / scale^-100 underflows
    println!("{:?}", w.inverse_cdf(1e-50)); // should be 0.0003088161777508183
    println!("{:?}", w.inverse_cdf(1e-20)); // should be 0.0006161692817189387
    // remedy: scale * (-ln(1 - p))^(1 / shape)
    println!("{:?}", 2f64.powi(-10) * (-(-1e-50f64).ln_1p()).powf(1.0 / 100.0));
}
```

This is the program that was run (release and debug mode give the same
output); it prints `0.0`, `inf`, `NaN`, `inf`, `0.0`,
`0.0006161727169554713` and, for the remedy, `0.0003088161777508183`.

## Expected

| shape, scale, p | Crate | mpmath (true value) | SciPy |
|---|---|---|---|
| 3.5, 2^-500, 0.5 | 0.0 | 2.7512099172225064e-151 | 2.7512099172225064e-151 |
| 3.5, 2^500, 0.5 | inf | 2.9479451064801336e150 | 2.9479451064801336e150 |
| 3.5, 2^500, 0 | NaN | 0 | 0 |
| 2.1, 2^500, 1e-5 | inf | 1.3615904136475406e148 | 1.3615904136475406e148 |
| 100, 2^-10, 1e-50 | 0.0 | 0.0003088161777508183 | 0.0003088161777508183 |
| 100, 2^-10, 1e-20 | 0.0006161727169554713 | 0.0006161692817189387 | 0.0006161692817189387 |

Second opinion: SciPy (`scipy.stats.weibull_min.ppf`) agrees with mpmath
to the last digit at every row; it passes all 1170 untagged points of the
dataset.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route, scale (-log1p(-p))^(1/shape) at
   60 extra bits (`_mp_weibull_ppf`), is tested against the cdf and sf
   (F(ppf(p)) = p) and against scale exp(log(-log(1 - p)) / shape) at 1200
   bits, in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: two random scans (Python `random.Random(20261017)`), each
   20,000 points: shape and scale as in finding 037 (shape 10^u, u uniform
   in [-1, 2.5]; scale 2^m), p half 10^u with u uniform in [-320, 0], a
   quarter uniform in (0, 1), a quarter 1 - 10^v with v uniform in
   [-16, 0], 6 significant digits; mpmath at 50 and 100 digits agreed at
   every point. With m in [-1000, 1000], `inverse_cdf` fails 11,970; with
   m in [-20, 20], 3,180. Every failure has scale^-shape or
   L / scale^-shape outside the normal range (checked point by point); on
   the grid no point with both inside fails.
5. **Root cause** (from the source, `src/distribution/weibull.rs`): `new`
   stores `scale_pow_shape_inv = scale.powf(-shape)`, and `inverse_cdf`
   computes `(-((-p).ln_1p() / scale_pow_shape_inv)).powf(1.0 / shape)`,
   so scale^shape is applied inside the power, where it can overflow or
   underflow although the answer is moderate. Confirmed with a bit-exact
   copy of the formula (every grid and random point gives the same bits as
   statrs): `scale * (-(-p).ln_1p()).powf(1.0 / shape)` passes every grid
   point (worst 2.8e-14). In the scans it fails only where
   (-ln(1 - p))^(1/shape) alone is below 2.2e-308 while the answer is not
   (206 and 6 points, all with scale above 1); scaling that power, for
   example through its exponent, would cover those too.

## Notes

Boost.Math's `weibull_distribution` computes the quantile as
`scale * pow(-boost::math::log1p(-p), 1 / shape)`; SciPy's `weibull_min`
likewise multiplies by the scale last.
