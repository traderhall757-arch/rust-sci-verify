# 024: statrs InverseGamma cdf and sf panic for x below 5.56e-309

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::InverseGamma`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** medium (a panic, which ends the calling thread, instead of
  the correct 0 or 1; no silently wrong value)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 30 of 531 points of
  `reference/data/invgamma_cdf.jsonl` and the same 30 of
  `reference/data/invgamma_sf.jsonl`: x = 5e-324 and x = 1e-310 at every one
  of the 15 shapes. The exact list is `INVGAMMA_TINY_X` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Every x from the
  smallest subnormal up to about 5.56e-309 (= 1 / f64::MAX) panics, at every
  shape (with rate 1; with another rate, every x below about
  rate * 5.56e-309). The true values there are 0 (`cdf`) and 1 (`sf`).

The other `InverseGamma::cdf` and `sf` failures (8 and 10 points, shapes
>= 1000 and two subnormal results) are finding 012's errors in `gamma_ur`
and `gamma_lr`, which these methods call; they are listed under
[finding 012](012-statrs-gamma-cdf.md). The quantile is
[finding 025](025-statrs-inverse-gamma-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, InverseGamma};

fn main() {
    let ig = InverseGamma::new(2.0, 1.0).unwrap();
    println!("{:?}", ig.cdf(5.6e-309)); // 0.0 (true value e^-1.8e308 rounds to 0)
    println!("{:?}", ig.sf(5.6e-309)); // 1.0
    println!("{:?}", 1.0 / 5.5e-309f64); // inf: 1/x overflows below 5.5626846462680035e-309
    println!("{:?}", ig.sf(5.5e-309)); // should be 1.0; panics
    println!("{:?}", ig.cdf(5.5e-309)); // should be 0.0; panics
}
```

This is the program that was run (release and debug mode): it prints the
first three lines and then panics with
``called `Result::unwrap()` on an `Err` value: XInvalid`` at
`src/function/gamma.rs:267`.

## Expected

| a | x | Method | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|---|
| 2 | 5.6e-309 | `cdf` | 0.0 | 0 | 0 |
| 2 | 5.6e-309 | `sf` | 1.0 | 1 | 1 |
| 2 | 5.5e-309 | `sf` | panic (XInvalid) | 1 | 1 |
| 2 | 5.5e-309 | `cdf` | panic (XInvalid) | 0 | 0 |
| every grid shape | 1e-310, 5e-324 | `cdf` | panic (XInvalid) | 0 | 0 |
| every grid shape | 1e-310, 5e-324 | `sf` | panic (XInvalid) | 1 | 1 |

Second opinion: SciPy (`scipy.stats.invgamma`) returns exactly 0 and 1 at
all 30 points of each dataset. The values need no precision: F(x) = Q(a, 1/x)
with 1/x > 1.8e308 is below e^-1e308, far below the smallest double.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route, Q(a, 1/x) and P(a, 1/x) with 1/x at
   60 extra bits (`_mp_invgamma_cdf`, `_mp_invgamma_sf`, on the gamma_sf and
   gamma_cdf series and continued fraction), is tested against
   `mpmath.gammainc` at 1/x, the closed forms for a = 1 (e^(-1/x)) and
   a = 1/2 (erfc(1/sqrt(x))), F + (1 - F) = 1, and, at x = 5e-324 and 1e-310,
   that the cdf is below 2^-1075 and the sf within 2^-54 of 1, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random points per method (Python
   `random.Random(20261015)`, shared with finding 025's scan): a = 10^u with
   u uniform in [-2, 4]; x = 10^v with v uniform in [-320, 300] (half) or
   x = 1 / (a + k sqrt(a)) with k uniform in [-30, 30] (half); all rounded to
   6 significant digits; mpmath at 50 and 100 digits (they agreed at all
   points). `cdf` panics at 210 points and `sf` at 213, every one with
   x <= 5.11e-309; no x >= 6.49e-309 panics. The other failures of the scan
   are finding 012's (section "Notes").
5. **Root cause** (from the source, `src/distribution/inverse_gamma.rs` and
   `src/function/gamma.rs`): `cdf` returns
   `gamma::gamma_ur(self.shape, self.rate / x)` and `sf` returns
   `gamma::gamma_lr(self.shape, self.rate / x)` for every finite x > 0. For
   x below rate / f64::MAX the quotient overflows to +inf; `checked_gamma_ur`
   and `checked_gamma_lr` reject x = +inf with `Err(XInvalid)`, and the
   unchecked `gamma_ur` / `gamma_lr` `unwrap` that error. Confirmed: guarding
   the quotient (returning 0 from `cdf` and 1 from `sf` when `rate / x` is
   infinite, which are the correctly rounded values there) removes all 423
   panics of the scan and changes no other result.

## Notes

Only rate 1 is checked. The off-grid scan's other
`cdf` / `sf` failures (532 and 533) are finding 012's: errors in normal
results for shapes from 643 (`cdf`) and 738 (`sf`) up, at most 1.8e-12 below
a = 1000 and 2.8e-11 overall (526 and 523 points), and subnormal results
returned as 0 (6 and 10 points; errors measured relative to 2.2e-308, up to
5.7e-3).
