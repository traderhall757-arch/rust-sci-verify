# 048: statrs Levy cdf and sf have about 10 correct digits, return 1 and 0 where x - location overflows, and lose the far upper tail of sf

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Levy`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location -2^1023 and
  scale 2^1023, `cdf(2^1023)` is 1 instead of 0.480 and `sf(2^1023)` is 0
  instead of 0.520; with location 0 and scale 2^-500, `sf(2^1023)` is 0
  instead of 4.7e-230. Elsewhere about 10 correct digits instead of 15 to 16)
- **Status:** approved (Claude, checklist complete, 2026-10-17)
- **Affected inputs:** 873 of 1776 points of `reference/data/levy_cdf.jsonl`
  and 381 of 1776 of `reference/data/levy_sf.jsonl` (tolerance 1e-13). The
  exact lists are `LEVY_CDF` and `LEVY_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). With
  t = (x - location) / scale and s = sqrt(1 / (2 t)), statrs computes
  `cdf` = erfc(s) and `sf` = erf(s). Four causes:
  1. **Finding 004's error** (869 `cdf` points, 371 `sf` points): statrs's
     erfc has about 10 correct digits for 0.5 <= s <= 26.5 (finding 004),
     and its erf is 1 minus the same erfc value for s >= 0.5 (one shared
     routine, `erf_impl`). `cdf` fails at every grid point with
     t = 7.2e-4 to 2 (s = 26.4 to 0.5; 9.1e-12 to 1.03e-10), on all six
     (location, scale) pairs. `sf` carries the same absolute error,
     diluted by a value near 1: it fails for t = 0.109 to 2 (s = 2.14 to
     0.5; 1.0e-13 to 9.5e-11).
  2. **Rounded s in the deep lower tail** (36 `cdf` points, all also in
     cause 1): erfc(s) has condition number about 2 s^2 in s, up to 1400
     where the result is still a normal double, so the two roundings in
     s = sqrt(0.5 c / (x - mu)) cost up to about 1.5e-13 there even with an
     exact erfc. A copy with only finding 004 fixed still fails these 36
     points (t = 7.6e-4 to 9.8e-4, results 5.3e-287 to 1.1e-224, 1.1e-13
     to 1.5e-13).
  3. **x - location above the largest double** (4 points each: location
     -2^1023, scale 2^1023, x = 2^1023 to 1.75 2^1023, so t = 2 to 2.75):
     `x - self.mu` overflows to inf although t is moderate, so s = 0 and
     `cdf` returns exactly 1, `sf` exactly 0 (finding 044's cause in
     `Laplace`, 046's in `Cauchy`).
  4. **0.5 c / (x - location) below the smallest normal double** (6 `sf`
     points, scale 2^-500, x = 3 2^550 to 2^1023, t = 3 2^1050 to 2^1523):
     the quotient is subnormal (digits lost: 6e-8 off at t = 3 2^1050, 6%
     at t = 3 2^1070) or 0 (from t = 2^1100, where `sf` returns 0 although
     the true value is 2.2e-166 to 4.7e-230). `cdf` is about 1 there and right.
  Every other point passes; the results 0 and 1 outside the support
  (x <= location) are right.

The quantile fails for another reason: [finding 049](049-statrs-levy-inverse-cdf.md).
Finding 004's error also reaches `Normal` ([finding 009](009-statrs-normal.md)).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Levy};

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let s = Levy::new(0.0, 1.0).unwrap(); // (location, scale)
    println!("{:?}", s.cdf(1.0)); // should be 0.3173105078629141
    println!("{:?}", s.cdf(0.0078125)); // should be 1.1224297172982926e-29
    println!("{:?}", s.sf(1.0)); // should be 0.6826894921370859
    let l = Levy::new(-big, big).unwrap();
    println!("{:?}", l.cdf(big)); // should be 0.4795001221869535 (t = 2)
    println!("{:?}", l.sf(big)); // should be 0.5204998778130465
    let n = Levy::new(0.0, 2f64.powi(-500)).unwrap();
    println!("{:?}", n.sf(3.0 * 2f64.powi(570))); // should be 4.0957340993776857e-162
    println!("{:?}", n.sf(big)); // should be 4.651557173616831e-230
    println!("{:?}", l.inverse_cdf(0.5)); // should be 1.0769164661541097e308 (finding 049)
    println!("{:?}", l.inverse_cdf(0.4999999999)); // should be 1.0769164652323077e308
}
```

This program was run in release and debug mode with the same output:
0.3173105078901144, 1.1224297172702965e-29, 0.6826894921098856, 1.0, 0.0,
4.3441820334104074e-162, 0.0, inf and inf.

## Expected

| Function (mu, c, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| cdf (0, 1, 1) | 0.3173105078901144 | 0.3173105078629141 | 0.31731050786291415 | 8.6e-11 |
| cdf (0, 1, 1/128) | 1.1224297172702965e-29 | 1.1224297172982926e-29 | 1.1224297172982928e-29 | 2.5e-11 |
| cdf (100, 8, 100 + 8 56/65536) | 1.750289210404745e-256 | 1.7502892104903754e-256 | 1.7502892104906315e-256 | 4.9e-11 |
| cdf (-2^1023, 2^1023, 2^1023) | 1.0 | 0.4795001221869535 | 1.0 | 1.1 |
| sf (0, 1, 1) | 0.6826894921098856 | 0.6826894921370859 | 0.6826894921370859 | 4.0e-11 |
| sf (0, 1, 1/4) | 0.9544997361056748 | 0.9544997361036416 | 0.9544997361036416 | 2.1e-12 |
| sf (-2^1023, 2^1023, 2^1023) | 0.0 | 0.5204998778130465 | 0.0 | 1 |
| sf (0, 2^-500, 3 2^550) | 4.1940319677465133e-159 | 4.19403171776275e-159 | 0.0 | 6.0e-8 |
| sf (0, 2^-500, 3 2^570) | 4.3441820334104074e-162 | 4.0957340993776857e-162 | 0.0 | 0.061 |
| sf (0, 2^-500, 2^1023) | 0.0 | 4.651557173616831e-230 | 0.0 | 1 |

Second opinion: SciPy (`scipy.stats.levy`) agrees with mpmath within a few
ulps at every cause-1 and cause-2 row, but is wrong at every cause-3 and
cause-4 row (it returns 1 and 0 where x - loc overflows, and 0 in the far
tail). The generator's route (erfc(s) and erf(s) with
s = sqrt(c / (2 (x - mu))) at 60 extra bits) is tested in
`reference/tests/test_functions.py` against quadrature of the density
(no error function; including t = 2 and 2.75 with location -2^1023, and the
deep lower tail at t = 46/65536), against the series of erf for the far
upper tail (t = 3 2^1050, 3 2^1070, 2^1200, 2^1523 with scale 2^-500, and
F = 1 - Q there), and outside the support.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the route is tested as described under
   "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: a scratch copy of statrs's Levy `cdf` and `sf` (with and
   without remedies) and mpmath at 50 and 100 digits (they agreed at all
   points). Scan 1 (Python `random.Random(20261028)`): 10,000 random
   (mu, c, x): c = 10^v with v uniform in [-300, 300], mu = 0 (a quarter)
   or +-10^u with u uniform in [-5, 5], both rounded to 6 significant
   digits, and x = mu + c 10^w rounded to a double, with w uniform in
   [-3.2, 1] for half, in [1, 300] for a quarter and in [-3.2, -3] (the
   deep lower tail) for a quarter (points with x <= mu or x infinite are
   redrawn). `cdf` fails 6,039, every one with s from 0.5 to 26.7
   (cause 1; above s = 26.5 the results are just below the smallest
   normal double, judged relative to 2.2e-308); `sf` fails 1,504, every
   one with s from 0.5 to about 2.2 (cause 1). A copy with only finding
   004's constants fixed still fails 407 `cdf` points, every one with
   s = 22.6 to 26.5, by 1.0e-13 to 1.7e-13 (cause 2), and no `sf` point.
   Scan 2 (`random.Random(20261029)`): 2,000 points with mu = -10^u,
   c = 10^v and x = 10^w, u, v and w uniform in [307.5, 308.25], 6
   significant digits. `cdf` and `sf` each fail 1,536: 845 where x - mu
   overflows (cause 3, exactly 1 and 0) and 691 through cause 1. Scan 3
   (`random.Random(20261030)`): 2,000 points with mu = 0, c = 10^v with v
   uniform in [-308, -200] and x = 10^w with w uniform in [100, 308]:
   `sf` fails 1,997, every one where 0.5 c / x is subnormal or 0
   (cause 4); `cdf` fails none.
5. **Root cause** (from the source, `src/distribution/levy.rs` and
   `src/function/erf.rs`): `cdf` is
   `erfc(((0.5 * self.c) / (x - self.mu)).sqrt())` and `sf` the same with
   `erf`. `erf_impl` returns erfc's value for s >= 0.5 and erf as
   1 minus it, so both carry finding 004's wrong interval constants
   (cause 1); erfc of a rounded square root magnifies the rounding
   (cause 2); `x - self.mu` overflows before the division brings it back
   (cause 3); and `0.5 * self.c / (x - self.mu)` underflows where its
   square root would still be a normal double (cause 4). Confirmed with a
   copy that uses finding 004's fix (each interval constant as
   `<constant>_f32 as f64`), forms s = sqrt(0.5 c) / sqrt(x - mu) (or from
   t = x / c - mu / c where x - mu overflows), and for `cdf` with s >= 0.5
   takes erfc's exp(-s^2) at u = 0.5 c / (x - mu) itself, corrected by
   u's rounding error (`(-u).mul_add(d, 0.5 * c) / d`), instead of at the
   rounded s^2. This copy passes every grid point (worst 4.4e-16 for both)
   and every point of the three scans (worst 6.8e-14 for `cdf`, in the
   deep lower tail where rounding u once costs up to 7.8e-14, and 4.2e-16
   for `sf`).

## Notes

The tolerance assumes a careful method rounds u = c / (2 (x - mu)) once:
F's condition number in u is about u, at most 698 at a normal result on
this grid, so that costs up to 7.8e-14 (tested in
`reference/tests/test_functions.py`). statrs's erfc(sqrt(u)) costs about
twice that through the square root (cause 2); on the grid that never
decides a point by itself, because cause 1 fails the same points first.
