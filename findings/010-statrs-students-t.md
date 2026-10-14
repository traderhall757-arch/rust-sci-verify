# 010: statrs StudentsT cdf and sf return 0 for huge |x|, 0.5 for tiny |x|, and lose digits for large df

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::StudentsT`,
  methods `cdf` and `sf` of `ContinuousCDF`; both call
  `statrs::function::beta::beta_reg`)
- **Severity:** high (a silently wrong answer: 0 instead of a normal result
  for |x| above 1.34e154 when df < about 1.85), plus an accuracy loss that
  grows with df: up to 8.9e-8 on the grid and 4.3e-4 off it (about 4 correct
  digits)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 144 of 1007 points of `reference/data/t_cdf.jsonl` and
  144 of 1007 points of `reference/data/t_sf.jsonl` (tolerance 2e-13); the
  exact lists are `T_CDF` and `T_SF` in `crates/check-statrs/src/expected.rs`,
  the report is [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  Only location 0 and scale 1 are checked. `sf` is the mirror image of
  `cdf`. Three groups, by cause (section "Root cause"):
  1. **0 for huge |x|** (3 points of each dataset): `cdf(-1e300)` for
     df = 0.1, 0.5 and 1 returns 0; the true values are 4.2e-31, 3.2e-151 and
     3.2e-301 (`sf(1e300)` likewise). Every |x| > 1.34e154 gives exactly 0 or 1.
     That is right only where the true lower tail is below the smallest
     double; on the random check below the largest df with a nonzero true
     value at such x is 1.85.
  2. **Exactly 0.5 for tiny |x|** (38 points): x = +-1e-10 gives 0.5 for every
     df on the grid; the error is |x| f(0) / 0.5, from 3.0e-11 (df = 0.1) to
     8.0e-11 (large df). Whenever x^2 / df < about 1e-15 the result is exactly
     0.5, and for x^2 / df up to about 1e-8 it is inaccurate; off the grid this
     reaches 4.3e-4 (df = 5.1e8, x = -5.4e-4 gives 0.5 for 0.49978).
  3. **Large df** (103 points, all with df >= 20): the error grows with df,
     worst per df 2.9e-13 (df = 20), 1.2e-12 (100), 1.8e-12 (1000),
     2.9e-11 (1e4), 5.1e-10 (1e5), 4.1e-9 (1e6) and 8.9e-8 (1e8), across the
     bulk and the tails (e.g. df = 1e5 fails at every x from -30 to 1.5).

`StudentsT::inverse_cdf` has separate, more severe problems:
[finding 011](011-statrs-students-t-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, StudentsT};

fn main() {
    let t = |df| StudentsT::new(0.0, 1.0, df).unwrap();
    println!("{:?}", t(1.0).cdf(-1e300)); // 0, should be 3.1830988618379065e-301
    println!("{:?}", t(0.1).cdf(-1e300)); // 0, should be 4.1738031371732015e-31
    println!("{:?}", t(1.0).sf(1e300)); // 0, should be 3.1830988618379065e-301
    println!("{:?}", t(5.0).cdf(1e-10)); // 0.5, should be 0.5000000000379606
    println!("{:?}", t(100.0).cdf(-1.5)); // should be 0.06838252906234443
    println!("{:?}", t(1e5).cdf(-3.0)); // should be 0.0013502304420323595
    println!("{:?}", t(1e8).cdf(-3.0)); // should be 0.0013498983640187472
    println!("{:?}", t(1e8).cdf(-0.01)); // should be 0.4960106436953425
    println!("{:?}", t(1e6).sf(2.0)); // should be 0.022750266925659603
}
```

This is the program that was run (release mode); it gives the "Crate"
column below.

## Expected

| df | x | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 1 | -1e300 | `cdf` | 0.0 | 3.1830988618379065e-301 | 0.0 | 1 |
| 0.1 | -1e300 | `cdf` | 0.0 | 4.1738031371732015e-31 | 0.0 | 1 |
| 1 | 1e300 | `sf` | 0.0 | 3.1830988618379065e-301 | 0.0 | 1 |
| 5 | 1e-10 | `cdf` | 0.5 | 0.5000000000379606 | 0.5000000000379607 | 7.6e-11 |
| 100 | -1.5 | `cdf` | 0.06838252906242925 | 0.06838252906234443 | 0.06838252906234442 | 1.2e-12 |
| 1e5 | -3 | `cdf` | 0.0013502304419257375 | 0.0013502304420323595 | 0.0013502304420323617 | 7.9e-11 |
| 1e8 | -3 | `cdf` | 0.0013498983815248369 | 0.0013498983640187472 | 0.0013498983640187468 | 1.3e-8 |
| 1e8 | -0.01 | `cdf` | 0.4960106877789744 | 0.4960106436953425 | 0.4960106436953424 | 8.9e-8 |
| 1e6 | 2 | `sf` | 0.022750266934293436 | 0.022750266925659603 | 0.02275026692565962 | 3.8e-10 |

Second opinion: SciPy (`scipy.stats.t`) agrees with mpmath within 1.7e-15 at
every row except the three at |x| = 1e300, where it returns 0 as well. There
the second route is the positive-term series DLMF 8.17.8 (`_mp_betainc`),
which gives the same doubles for df = 0.1, 0.5 and 1, and for df = 1 the
closed form F(x) = atan(-1/x) / pi (tested in
`reference/tests/test_functions.py`).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route F(-|x|) = I_h(df/2, 1/2) / 2 by the
   continued fraction DLMF 8.17.22 (`_mp_betainc_cf`) is tested against the
   series DLMF 8.17.8 and `mpmath.betainc`, against closed forms for
   df = 1, 2 and 3, against the hypergeometric form
   1/2 + x c 2F1(1/2, (df+1)/2; 3/2; -x^2/df) (at 300 digits, including
   df = 1e8, x = -30), and cdf + sf = 1, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random points (Python `random.Random(20261014)`):
   df = 10^u with u uniform in [-1.3, 9], x = +-10^v with v uniform in
   [-12, 305], each rounded to 6 significant digits; mpmath at 50 and
   100 digits (they agreed at all 20,000). `cdf` fails 1,357 (6.8%). By the
   fixes of the experiment below: making huge |x| work leaves 496, then
   exact 1 - h and no cut-off leave 135, then `ln_1p` leaves 82, then an
   accurate prefactor leaves 8 (worst 6.9e-9 at df near 1e8, the continued
   fraction's own rounding). So every failure has one of the causes below.
   Worst error per decade of df: 4.6e-8 (10 to 100), 3.6e-7 (100 to 1000),
   4.2e-7, 3.8e-6, 5.9e-6, 1.7e-5, 1.6e-4 and 4.3e-4 (1e8 to 1e9); the
   largest are at small |x| (cause 2). All 681 random points with
   x < -1.34e154 and a nonzero true value (df up to 1.85) give 0.
5. **Root cause** (confirmed by experiment): `StudentsT::cdf` computes
   `h = df / (df + k * k)` and `0.5 * beta_reg(df / 2, 0.5, h)`
   (`src/distribution/students_t.rs`), and `beta_reg`
   (`src/function/beta.rs`) is the Numerical Recipes continued fraction with
   the prefactor `exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a ln x +
   b ln(1 - x))`. A copy of both functions (scratch, not committed; it
   reproduces statrs bit for bit on all 1007 grid and 20,000 random points)
   with switchable fixes shows:
   - **Cause 1:** `k * k` overflows for |k| > 1.34e154, so h = 0 and the
     result is 0 (or 1). Computing the prefactor in log space,
     ln h = ln df - 2 ln|x| - ln_1p(df/x^2), fixes all huge-|x| points
     (worst then 8.9e-15).
   - **Cause 2:** `beta_reg` gets only h and forms 1 - h itself. For tiny
     x^2 / df, h rounds to 1 or close to it, so 1 - h has lost its digits,
     and `beta_reg` sets the prefactor to 0 when |h - 1| <= 1e-15 (`ulps_eq!`
     with `MODULE_EPS`), which gives exactly 0.5. Passing
     z = x^2 / (df + x^2) separately and removing the cut-off fixes all 38
     points (df <= 10 then within 5.1e-14 everywhere on the grid).
   - **Cause 3:** for large df three losses add up. (a) The prefactor
     subtracts `ln_gamma` values of size about (df/2) ln(df/2) (8.4e8 at
     df = 1e8), each with an absolute error of up to about 1e-13 (measured:
     statrs `ln_gamma(25)` is off by 1.1e-13), and the result is
     exponentiated. Replacing it by the correctly rounded ln B(df/2, 1/2)
     removes the df = 20 to 100 failures (worst then 3.3e-14) and cuts the
     worst error at df = 1e4 from 1.3e-11 to 6.9e-13 and at df = 1e8 from
     5.5e-8 to 1.2e-8 (after the cause 2 fix). (b) `a * ln(h)` with
     a = df/2 amplifies the rounding of h by df/2; using
     ln h = -ln_1p(x^2/df) instead cuts the df = 1e5 error further from
     1.3e-11 to 2.1e-12 and the df = 1e8 error from 1.2e-8 to 2.9e-9. (c) The continued fraction itself, evaluated in double
     precision with h near 1, loses the rest: run with exact arithmetic on
     the same rounded inputs it is good to 4e-13 at df = 1e5 and 1.2e-10 at
     df = 1e8, but in double precision it gives exactly the remaining
     2.1e-12 and 2.9e-9 (mpmath 53-bit emulation matches the Rust copy).
     With all of these fixes the grid has 25 failures left, all at
     df >= 1e4 (worst 2.9e-9 at df = 1e8), which only a different method
     for large a (Boost uses series and asymptotic expansions there) removes.

## Notes

Boost.Math's `students_t` cdf also uses I_h(df/2, 1/2) but passes both h and
1 - h to `ibeta`/`ibetac` (the smaller one computed directly), computes the
prefactor with a Lanczos approximation that combines the gamma terms before
rounding, and has special methods for large a with small b. For huge |x|
SciPy is also wrong, so the report lists those points as "both wrong".
