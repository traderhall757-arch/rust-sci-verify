# 025: statrs InverseGamma inverse_cdf off by up to 46% for shapes near 108 to 170, and in the far tails

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::InverseGamma`,
  method `inverse_cdf` of `ContinuousCDF`; the search is
  `distribution::internal::newton_raphson_quantile`)
- **Severity:** high (silently wrong answers: `InverseGamma(150, 1)` gives
  the median as 0.005859375 instead of 0.0066815, 12% off, and 42% off at
  p = 1e-10)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 38 of 405 points of
  `reference/data/invgamma_ppf.jsonl` (tolerance 2e-13). The exact list is
  `INVGAMMA_PPF` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Only rate 1 is
  checked. Three groups, by cause (section "Root cause"):
  1. **Shapes from about 108 to 170** (20 grid points, all at a = 150, every
     p from 1e-100 to 0.999; errors 1.8e-4 to 42%): the density overflows,
     and the search stops at its first guess, a power of two or the midpoint
     of two (0.005859375 = 3/512 for every p from 1e-5 to 0.95). Also 1 point
     at a = 1e4, p = 1/2 + 1e-10 (3.1e-13), where the density is NaN.
  2. **Upper tail for small shapes** (2 grid points: a = 0.1 and 0.5 at
     p = 1 - 1e-15; 1.4% and 0.6% off; off the grid a up to 1.8 and
     p >= 1 - 5.9e-14, up to 77% off): the bracket is built from the
     rounded `cdf`, misses the root, and the search returns its end.
  3. **Subnormal p** (15 grid points: p = 5e-324 at every shape; 0.6% to
     3.8% off): `cdf` near p has one significant bit.

Group 1 is the serious one: for a in that range the quantile is wrong at
almost every p. The same search is used by `Chi::inverse_cdf` (not yet
checked).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Continuous, ContinuousCDF, InverseGamma};

fn main() {
    let ig = |a: f64| InverseGamma::new(a, 1.0).unwrap();
    println!("{:?}", ig(150.0).inverse_cdf(0.5)); // should be 0.006681508583437417
    println!("{:?}", ig(150.0).inverse_cdf(0.01)); // should be 0.005556999974976799
    println!("{:?}", ig(150.0).inverse_cdf(1e-100)); // should be 0.0017296089000173908
    println!("{:?}", ig(150.0).pdf(0.005)); // should be 1.0371243066031393
    println!("{:?}", ig(0.1).inverse_cdf(0.999999999999999)); // should be 1.6598386301668927e150
    println!("{:?}", ig(0.5).inverse_cdf(0.999999999999999)); // should be 1.2752773318461996e30
    println!("{:?}", ig(1.0).inverse_cdf(5e-324)); // should be 0.001343291471963653
    println!("{:?}", ig(1e4).inverse_cdf(0.5000000001)); // should be 0.0001000033334249437
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| a | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| 150 | 0.5 | 0.005859375 | 0.006681508583437417 | 0.0066815085834374165 | 0.12 |
| 150 | 0.01 | 0.005859375 | 0.005556999974976799 | 0.005556999974976799 | 5.4e-2 |
| 150 | 1e-100 | 0.00146484375 | 0.0017296089000173908 | 0.0017296089000173906 | 0.15 |
| 150 | (`pdf` at x = 0.005) | inf | 1.0371243066031393 | | |
| 0.1 | 1 - 1e-15 | 1.636695303948071e150 | 1.6598386301668927e150 | 1.659838630166961e150 | 1.4e-2 |
| 0.5 | 1 - 1e-15 | 1.2676506002282294e30 | 1.2752773318461996e30 | 1.275277331846191e30 | 6.0e-3 |
| 1 | 5e-324 | 0.0013926154463651573 | 0.001343291471963653 | 0.0013432914719636532 | 3.7e-2 |
| 1e4 | 1/2 + 1e-10 | 0.00010000333342491302 | 0.0001000033334249437 | 0.0001000033334249437 | 3.1e-13 |

(1 - 1e-15 is the double 0.999999999999999; 1.636695303948071e150 is
2^499 and 1.2676506002282294e30 is 2^100.) The density's true value is
x^(-a-1) e^(-1/x) / Gamma(a) evaluated in mpmath through its logarithm.

Second opinion: SciPy (`scipy.stats.invgamma.ppf`) agrees with mpmath within
4.1e-14 at every row. Off these rows it fails 11 of the 15 p = 5e-324 points
(all but a = 1, 1.5, 2 and 3); there the second route is F(ppf(p)) = p
checked with `mpmath.gammainc`, and for a = 1 the closed form: F(x) =
e^(-1/x), so the quantile is -1 / log(p), which gives the mpmath value.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route (`_mp_invgamma_ppf`) is 1/t with t
   the gamma quantile on the upper tail Q(a, t) = p (p <= 1/2) or the lower
   tail P(a, t) = 1 - p (Newton on the logarithm, shared with `gamma_ppf`).
   It is tested against F(ppf(p)) = p by `mpmath.gammainc` (a different
   method) at the points of every group above, against the closed form
   -1 / log(p) for a = 1 from p = 5e-324 to 1 - 1e-15, and ppf(0) = 0, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261015)`, after
   finding 024's 40,000): a = 10^u with u uniform in [-2, 4]; p = 10^v with v
   uniform in [-323, log10(0.5)] (40%), uniform in (0, 1) (30%) or
   1 - 10^v with v uniform in [-16, log10(0.5)] (30%); all rounded to 6
   significant digits; mpmath at 50 and 100 digits (they agreed at all
   points). 54 quantiles above the largest double are left out. Of the
   2,946 others, 97 fail (3.3%):
   - 49 with a from 108.0 to 170.0 (half of the 99 points in that range),
     up to 46% off (a = 164.9, p = 0.99971: 0.01171875 instead of 0.0080482);
   - 11 with a from 0.045 to 1.82 and p >= 1 - 5.9e-14, up to 77% off
     (a = 0.0903, p = 0.9999999999999997: 9.66e170 instead of 4.12e171);
   - 35 with subnormal p, up to 3.7% off, at every range of a;
   - 2 just above the tolerance: 2.4e-13 (a = 9886.25, p = 0.280135) and
     2.1e-13 (a = 0.0108, p = 0.0268).
5. **Root cause** (confirmed with a copy of `newton_raphson_quantile`, bit
   for bit equal to statrs on all 2,946 scan points, with switchable fixes):
   - **Cause 1, the density overflows.** `InverseGamma::pdf` computes
     `rate.powf(shape) * x.powf(-shape - 1.0) * (-rate / x).exp() / gamma(shape)`.
     The search evaluates it at x = 2^-k well below the quantile, where
     x^(-a-1) exceeds 1.8e308 once a is above about 108, so the density is
     +inf (or NaN, inf times an underflowed exponential, and for every
     a > 171.6, where `gamma(a)` overflows). With an infinite density the
     Newton step is x - residual / inf = x, which passes the test
     `(newton - x).abs() <= accuracy * x.abs()`, so the search returns its
     current iterate, the first midpoint of the bracket. A NaN density only
     forces bisection on the rounded `cdf`, which costs accuracy at
     a = 1e4 (finding 012's error in `cdf`). Forming the density from its
     logarithm, `exp(-(a + 1) ln x - 1/x - ln_gamma(a))`, removes all 21 grid
     points of this group and 50 of the scan's failures.
   - **Cause 2, the bracket.** The search first doubles `high` while
     `cdf(high) < p`. In the upper tail `cdf` is 1 - sf rounded to a double,
     with an absolute error of up to 1.1e-16, so for p = 1 - 1e-15 it can
     already reach p at a point below the quantile: for a = 0.1 the
     bracket stops at [2^498, 2^499] while the quantile is 1.66e150 > 2^499.
     The Newton phase keeps the bracket as an invariant and converges to its
     end. The search itself solves sf(x) = 1 - p in the upper half; building
     the bracket from `sf` there too removes both grid points of this group
     and 11 of the scan's failures, and changes no value that was within
     tolerance.
   - **Cause 3, subnormal p.** At p = 5e-324 the root lies where the true
     cdf is subnormal, so the computed `cdf` takes only a few values
     (multiples of 5e-324, or 0 where finding 012 flushes it) over a range
     of x about 5e-4 wide in relative terms, and the residual `cdf(x) - p`
     cannot locate the root. Even with exact closed forms for a = 1
     (`cdf` = e^(-1/x), `sf` = -expm1(-1/x)) the same search returns
     1.3433522593760319e-3 instead of 1.343291471963653e-3 (4.5e-5 off).
     Only a search on log F avoids this; it is finding 013's cause 4 for the
     gamma quantile. With causes 1 and 2 fixed, the 15 grid points and 35
     scan points of this group are all that fail, plus the 2 near-tolerance
     scan points above (2.0e-13 and 2.1e-13).

## Notes

`InverseGamma::pdf` itself is public, so the overflow in cause 1 also gives
direct callers +inf or NaN for shapes above about 108 (seen in the
reproduction; the density is not checked against reference data in this
project). The bracket and the stopping rule are in the shared
`newton_raphson_quantile`, so `Chi::inverse_cdf`, its only other caller,
should be checked next. Boost.Math computes this quantile as
1 / gamma_q_inv(a, p) with no density of its own; SciPy as
1 / gammainccinv(a, p).
