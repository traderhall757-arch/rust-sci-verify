# 023: statrs Exp cdf returns 0 for x below 5.6e-17 and loses digits for x below about 5e-4

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Exp`, method
  `cdf` of `ContinuousCDF`)
- **Severity:** high (a silently wrong answer: `cdf(1e-20)` is 0 instead of
  1e-20, and every x from 1e-310 to about 5.6e-17 gives 0), plus an
  accuracy loss of about 1.1e-16 / x for larger x: 3.1e-4 at x = 1e-13,
  8.3e-8 at 1e-10, 5.7e-13 at 1e-5
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 45 of 325 points of `reference/data/expon_cdf.jsonl`
  (tolerance 1e-13): every grid x from 1e-310 to 1e-5 (the decades
  1e-300 .. 1e-30, 10^-k for k = 5..20, and 1e-310). The exact list is
  `EXP_CDF` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Only rate 1 is
  checked; for another rate the same happens at rate * x.
  - **x below 2^-54 = 5.55e-17**: `cdf` returns exactly 0 (33 grid points
    from 1e-310 to 1e-17). At x = 5e-324 the true value is the smallest
    subnormal, and 0 is within tolerance under the subnormal rule
    (`docs/conventions.md`), so that point passes.
  - **x from 5.55e-17 to about 5e-4**: the result has an absolute error of
    up to about 1.1e-16, so a relative error of about 1.1e-16 / x: 11% at
    x = 1e-16, 3.1e-4 at 1e-13, 8.3e-8 at 1e-10, 5.7e-13 at 1e-5 (12 grid
    points). Off the grid, failures reach x = 4.5e-4 (1.2e-13).

`Exp::sf` (exact at every grid point) and `Exp::inverse_cdf` (worst 2.1e-16)
are correct.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Exp};

fn main() {
    let exp = Exp::new(1.0).unwrap();
    println!("{:?}", exp.cdf(1e-20)); // should be 1e-20
    println!("{:?}", exp.cdf(5e-17)); // should be 5e-17
    println!("{:?}", exp.cdf(1e-16)); // should be 1e-16
    println!("{:?}", exp.cdf(1e-13)); // should be 9.9999999999995e-14
    println!("{:?}", exp.cdf(1e-10)); // should be 9.999999999500001e-11
    println!("{:?}", exp.cdf(1e-5)); // should be 9.999950000166666e-6
    println!("{:?}", exp.cdf(0.000452551)); // should be 0.0004524486142417077
    println!("{:?}", -(-1e-10f64).exp_m1()); // -expm1(-x): 9.999999999500001e-11
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below. The last line shows the
remedy.

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 1e-310 | 0.0 | 1e-310 | 1.00000000000005e-310 | 4.5e-3 (relative to 2.2e-308) |
| 1e-20 | 0.0 | 1e-20 | 1e-20 | 1 |
| 5e-17 | 0.0 | 5e-17 | 5e-17 | 1 |
| 1e-16 | 1.1102230246251565e-16 | 1e-16 | 1e-16 | 0.11 |
| 1e-13 | 1.000310945187266e-13 | 9.9999999999995e-14 | 9.9999999999995e-14 | 3.1e-4 |
| 1e-10 | 1.000000082740371e-10 | 9.999999999500001e-11 | 9.999999999500001e-11 | 8.3e-8 |
| 1e-5 | 9.999950000172397e-6 | 9.999950000166666e-6 | 9.999950000166666e-6 | 5.7e-13 |
| 0.000452551 | 0.00045244861424176275 | 0.0004524486142417077 | 0.0004524486142417077 | 1.2e-13 |

Second opinion: SciPy (`scipy.stats.expon.cdf`) agrees with mpmath to the
last digit at every row except 1e-310, where it is one subnormal step off
(2.2e-16 relative to 2.2e-308). For small x the series
x - x^2/2 + x^3/6 - ... gives the mpmath value too (tested in
`reference/tests/test_functions.py`).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route, `-mpmath.expm1(-x)`
   (`_mp_expon_cdf`), is tested against the gamma_cdf route with shape 1
   (series DLMF 8.7.1), `mpmath.gammainc(1, 0, x)`, the Taylor series for
   small x, and cdf + sf = 1, in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random x (Python `random.Random(20261015)`): half
   x = 10^u with u uniform in [-320, 3], half uniform in (0, 750); each
   rounded to 6 significant digits; mpmath at 50 and 100 digits (they agreed
   at all points). `cdf` fails 9,805 (49%), every one with x <= 4.53e-4;
   9,433 of them return 0 (all with x < 5.55e-17). Worst error per decade of
   x: 1 up to 1e-17, 0.29 (1e-16), 3.3e-3 (1e-14), 2.8e-5 (1e-12), 4.6e-7
   (1e-10), 4.3e-9 (1e-8), 3.7e-11 (1e-6), 5.1e-13 (1e-4). `sf` fails none
   (worst 2.2e-16). `inverse_cdf` fails none of 20,000 random p (half
   10^u with u in [-320, 0], a quarter uniform in (0, 1), a quarter
   1 - 10^v with v in [-16, 0]; worst 2.2e-16).
5. **Root cause** (from the source, `src/distribution/exponential.rs`):
   `cdf` computes `1.0 - (-self.rate * x).exp()`. For small x, e^-x is close
   to 1 and is rounded to a double with an absolute error of up to 2^-54
   (half an ulp below 1); the subtraction is exact but keeps that absolute
   error, so the relative error is about 2^-54 / x, and for x < 2^-54,
   e^-x rounds to exactly 1 and the result is 0. Confirmed:
   `-(-self.rate * x).exp_m1()` (keeping the `x < 0` branch), the
   counterpart of the `ln_1p` that statrs already uses in `inverse_cdf`,
   passes every grid point (325 of 325, worst 1.3e-16) and all 20,000
   random points (worst 1.7e-16).

## Notes

`sf` returns `(-rate * x).exp()`, the correctly computed small tail, so
users who need 1 - F can use it. Boost.Math's `exponential_distribution`
computes the cdf as `-boost::math::expm1(-x * lambda)`; SciPy's
`expon.cdf` as `-expm1(-x)`.
