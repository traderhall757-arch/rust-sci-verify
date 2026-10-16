# 037: statrs Weibull cdf and sf give NaN, 0 or 1 where scale^-shape or x^shape leaves the double range

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Weibull`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with shape 3.5 and scale
  2^-500, `cdf(scale)` is NaN instead of 1 - 1/e = 0.632, and with shape 100
  and scale 1/1024, `cdf(9.765625e-5)` is 0 instead of 1e-100)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 177 of 1139 points of `reference/data/weibull_cdf.jsonl`
  and 180 of 1139 of `reference/data/weibull_sf.jsonl` (tolerance 1e-13).
  The exact lists are `WEIBULL_CDF` and `WEIBULL_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Every failure
  has an intermediate outside the normal double range (2.2e-308 to
  1.8e308), although t = (x / scale)^shape and the answer are not:
  1. **scale^-shape overflows to inf** (cdf 88, sf 94 points; scale 2^-500
     with shape >= 2.1): NaN where x = 0 or x^shape underflows (0 * inf),
     otherwise 1 (cdf) and 0 (sf) at every x > 0.
  2. **scale^-shape underflows to 0** (cdf 72, sf 71; scale 2^500 with
     shape >= 3.5): NaN where x^shape overflows (0 * inf), otherwise 0
     (cdf) and 1 (sf) at every x.
  3. **scale^-shape is subnormal** (cdf 9, sf 15; scale 2^500, shape 2.1:
     2^-1050): digits lost, or 1 / 0 where x^shape overflows.
  4. **x^shape underflows to 0 or a subnormal** (cdf only, 8 points; shape
     100 with scale 2^-10 and shape 1.5 or 2 with scale 2^-500): 0
     instead of a tiny cdf (1e-100), or up to 5.6e-4 off.

  Ordinary-looking parameters are affected too: with scale 2^-10 and shape
  100, `cdf` is wrong for every x below 7.8e-4, where the true value is
  below about 1e-10 (a scan of 400 x from 3.2e-4 to 2e-3). In the random scan
  below with scales from 2^-20 to 2^20, 21% of cdf and 13% of sf points fail.
  One further sf point in that scan (not on the grid) fails by rounding
  alone: computing x^shape and scale^-shape apart adds two roundings to t,
  which the sf's condition number t amplifies (1.2e-13 at t = 611).

`Weibull::inverse_cdf` fails for the same reason:
[finding 038](038-statrs-weibull-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Weibull};

fn main() {
    let tiny = 2f64.powi(-500); // 3.054936363499605e-151
    let huge = 2f64.powi(500); // 3.273390607896142e150
    let w = Weibull::new(3.5, tiny).unwrap(); // scale^-3.5 = 2^1750 overflows
    println!("{:?}", w.cdf(tiny)); // should be 0.6321205588285577
    println!("{:?}", w.cdf(0.0)); // should be 0.0
    println!("{:?}", w.sf(tiny)); // should be 0.36787944117144233
    let w = Weibull::new(3.5, huge).unwrap(); // scale^-3.5 = 2^-1750 underflows to 0
    println!("{:?}", w.cdf(huge)); // should be 0.6321205588285577
    let w = Weibull::new(2.1, huge).unwrap(); // scale^-2.1 = 2^-1050 is subnormal
    println!("{:?}", w.cdf(huge)); // should be 0.6321205588285577
    let w = Weibull::new(2.0, tiny).unwrap(); // x^2 is subnormal
    println!("{:?}", w.cdf(3.0549363634996048e-161)); // should be 1.0000000000000001e-20
    let w = Weibull::new(2.1, tiny).unwrap();
    println!("{:?}", w.sf(6.955155764017214e-150)); // should be 2.0061323053314493e-308
    let w = Weibull::new(100.0, 2f64.powi(-10)).unwrap(); // x^100 underflows to 0
    println!("{:?}", w.cdf(9.765625e-5)); // should be 1.0000000000000056e-100
    // remedy: t = (x / scale)^shape first (x / scale = 0.1 exactly)
    println!("{:?}", -(-(9.765625e-5f64 / 2f64.powi(-10)).powf(100.0)).exp_m1());
}
```

This is the program that was run (release and debug mode give the same
output); it prints `NaN` (four times), `1.0`, `1.0005576689441423e-20`,
`0.0`, `0.0` and, for the remedy, `1.0000000000000056e-100`.

## Expected

| shape, scale, x | Function | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|
| 3.5, 2^-500, 2^-500 | cdf | NaN | 0.6321205588285577 | 0.6321205588285577 |
| 3.5, 2^-500, 0 | cdf | NaN | 0 | 0 |
| 3.5, 2^-500, 2^-500 | sf | NaN | 0.36787944117144233 | 0.36787944117144233 |
| 3.5, 2^500, 2^500 | cdf | NaN | 0.6321205588285577 | 0.6321205588285577 |
| 2.1, 2^500, 2^500 | cdf | 1.0 | 0.6321205588285577 | 0.6321205588285577 |
| 2, 2^-500, 3.0549363634996048e-161 | cdf | 1.0005576689441423e-20 | 1.0000000000000001e-20 | 1.0000000000000001e-20 |
| 2.1, 2^-500, 6.955155764017214e-150 | sf | 0.0 | 2.0061323053314493e-308 | 2.0061323053315337e-308 |
| 100, 2^-10, 9.765625e-5 | cdf | 0.0 | 1.0000000000000056e-100 | 1.0000000000000056e-100 |

Second opinion: SciPy (`scipy.stats.weibull_min`) agrees with mpmath to
the last digit at every row except the subnormal sf (3.8e-14 relative to
2.2e-308, within tolerance); it passes all 1139 points of both datasets.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The routes, `-expm1(-t)` and `exp(-t)` with
   t = (x / scale)^shape at 60 extra bits (`_mp_weibull_cdf`,
   `_mp_weibull_sf`), are tested against the exponential distribution
   (shape 1), quadrature of the density, and cdf + sf = 1, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: two random scans (Python `random.Random(20261017)`), each
   20,000 points per function: shape 10^u with u uniform in [-1, 2.5]
   (4 significant digits), scale 2^m with m a uniform integer, x = v scale
   with v = 10^(w / shape), w uniform in [-300, 3], 6 significant digits
   (so t runs from 1e-300 to 1e3 and x / scale is exact); mpmath at 50 and
   100 digits agreed at every point. With m in [-1000, 1000], `cdf` fails
   15,274 and `sf` 14,140; with m in [-20, 20], 4,237 and 2,530. In both
   scans every failure but the one sf point above has scale^-shape or
   x^shape outside the normal range (checked point by point); on the grid no
   point with both inside fails.
5. **Root cause** (from the source, `src/distribution/weibull.rs`): `new`
   stores `scale_pow_shape_inv = scale.powf(-shape)`, and `cdf` and `sf`
   compute `-(-x.powf(shape) * scale_pow_shape_inv).exp_m1()` and
   `(-x.powf(shape) * scale_pow_shape_inv).exp()`. Each power can leave the
   double range while their product t = (x / scale)^shape is moderate.
   Confirmed with a bit-exact copy of both formulas (every grid and random
   point gives the same bits as statrs): computing t as
   `(x / scale).powf(shape)` instead passes every grid point (worst 2.2e-16
   cdf, 5.4e-14 sf) and every random point of both scans (worst 2.2e-16,
   5.3e-14).

## Notes

Boost.Math's `weibull_distribution` computes the cdf as
`-expm1(-pow(x / scale, shape))`, and SciPy's `weibull_min` standardises
x / scale first too. `pdf` and `ln_pdf` use the same stored power and very
likely fail the same way; they are not checked here.
