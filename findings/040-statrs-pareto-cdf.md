# 040: statrs Pareto cdf loses digits just above the scale, and cdf and sf go wrong where scale / x underflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Pareto`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with shape 2^-5 and scale 1,
  `cdf(1.0000000000000002)` is 0 instead of 6.9e-18; with shape 2^-5 and
  scale 2^-500, `sf(4.15e180)` is 0 instead of 4.5e-11), plus an accuracy
  loss of about 1.1e-16 / F(x) just above the scale
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 121 of 1062 points of `reference/data/pareto_cdf.jsonl`
  and 9 of 1062 of `reference/data/pareto_sf.jsonl` (tolerance 1e-13). The
  exact lists are `PARETO_CDF` and `PARETO_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes:
  1. **x just above the scale x_m** (`cdf` only, 120 grid points, x / x_m
     from 1 + 2^-52 to 1.001, every shape): 1 - (x_m/x)^alpha subtracts two
     nearly equal numbers, so the absolute error of about 1.1e-16 becomes a
     relative error of about 1.1e-16 / F(x). Where x_m / x rounds to 1 the
     result is exactly 0 (15 grid points). Ordinary parameters are hit:
     with shape 1 and scale 1, `cdf(1.0000000001)` has 9 correct digits.
  2. **x_m / x below 2.2e-308** (`sf` 9 points, `cdf` 1; scale 2^-500 with
     x / x_m = 1.1 2^1060, 2^1100 or 2^1500 and shapes 2^-5, 0.1, 0.5): the
     quotient is subnormal (digits lost) or 0, although (x_m/x)^alpha is a
     moderate number for a shape below 1. `sf` returns 0 or a value with 5
     correct digits; `cdf` returns 1 where the true value is
     1 - 4.5e-11.

`Pareto::inverse_cdf` fails elsewhere: [finding 041](041-statrs-pareto-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Pareto};

fn main() {
    let tiny = 2f64.powi(-500); // 3.054936363499605e-151
    // statrs takes (scale, shape).
    let p = Pareto::new(1.0, 2f64.powi(-5)).unwrap();
    println!("{:?}", p.cdf(1.0000000000000002)); // should be 6.9388939039072276e-18
    let p = Pareto::new(1.0, 1.0).unwrap();
    println!("{:?}", p.cdf(1.0000000001)); // should be 1.000000082640371e-10
    let p = Pareto::new(1.0, 100.0).unwrap();
    println!("{:?}", p.cdf(1.00001)); // should be 0.000999495171662341
    let p = Pareto::new(tiny, 2f64.powi(-5)).unwrap();
    println!("{:?}", p.cdf(4.149515568880993e180)); // should be 0.9999999999551158
    println!("{:?}", p.sf(4.149515568880993e180)); // should be 4.488424239121205e-11
    println!("{:?}", p.inverse_cdf(0.9999999999)); // should be 3.054928274980468e169
    let p = Pareto::new(tiny, 0.1).unwrap();
    println!("{:?}", p.sf(4.151358667303696e168)); // should be 1.2209030849522242e-32
    let p = Pareto::new(tiny, 0.5).unwrap();
    println!("{:?}", p.sf(1.0715086071862673e301)); // should be 1.688508503057271e-226
    let p = Pareto::new(1e-20, 0.5).unwrap();
    println!("{:?}", p.sf(1e300)); // should be 1e-160
    // Remedies: log1p of the exact x - x_m, and half the power applied twice.
    let (x, xm, a) = (1.0000000000000002f64, 1.0f64, 2f64.powi(-5));
    println!("{:?}", -(-a * ((x - xm) / xm).ln_1p()).exp_m1()); // 6.9388939039072276e-18
    let h = (1.0f64 - 0.9999999999).powf(-0.5 / a);
    println!("{:?}", tiny * h * h); // 3.054928274980469e169: one ulp from the true value
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below and the `inverse_cdf` row of
finding 041.

## Expected

| Function (alpha, x_m, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| cdf (2^-5, 1, 1.0000000000000002) | 0.0 | 6.9388939039072276e-18 | 0.0 | 1 |
| cdf (1, 1, 1.0000000001) | 1.000000082740371e-10 | 1.000000082640371e-10 | 1.000000082740371e-10 | 1.0e-9 |
| cdf (100, 1, 1.00001) | 0.0009994951716614242 | 0.000999495171662341 | 0.0009994951716623124 | 9.2e-13 |
| cdf (2^-5, 2^-500, 4.149515568880993e180) | 1.0 | 0.9999999999551158 | 1.0 | 4.5e-11 |
| sf (2^-5, 2^-500, 4.149515568880993e180) | 0.0 | 4.488424239121205e-11 | 0.0 | 1 |
| sf (0.1, 2^-500, 4.151358667303696e168) | 1.2209068108015856e-32 | 1.2209030849522242e-32 | 0.0 | 3.1e-6 |
| sf (0.5, 2^-500, 1.0715086071862673e301) | 0.0 | 1.688508503057271e-226 | 0.0 | 1 |
| sf (0.5, 1e-20, 1e300), off the grid | 9.99994433575849e-161 | 1e-160 | 0.0 | 5.6e-6 |

Second opinion: SciPy (`scipy.stats.pareto`) has the same two weaknesses
(it is right only in the third row), so every row was also computed by a
second mpmath route: (x_m/x)^alpha and 1 - (x_m/x)^alpha directly at 1200
bits, which agree with the table to every printed digit. In
`reference/tests/test_functions.py` the generator's route is also tested
against the exponential distribution (log(X / x_m) is exponential with rate
alpha, logarithm at 1200 bits), quadrature of the density, and the exact
value 2^(-1500 alpha) at x / x_m = 2^1500.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route (`_mp_pareto_cdf` and
   `_mp_pareto_sf`: -expm1(-t) and exp(-t) with t = alpha log(x / x_m)) is
   tested as described under "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (alpha, x_m, x) (Python
   `random.Random(20261019)`): alpha = 10^u with u uniform in [-1.6, 2] and
   x_m = 10^v with v uniform in [-300, 300], both rounded to 6 significant
   digits (so most scales are not powers of two); half with
   x = x_m (1 + 10^w), w uniform in [-16, 0], half with x = x_m 10^w, w
   uniform in [0, 620], kept where x is a normal double; mpmath at 50 and
   100 digits (they agreed at all points). `cdf` fails 7,937 (40%): 7,869
   with x / x_m below 1.019 (cause 1; 155 of them return 0) and 68 with
   x_m / x below 2.2e-308 (cause 2); worst error 1.4. `sf` fails 989, every
   one with x_m / x below 2.2e-308 (shapes up to 0.94, true values up to
   5.4e-9).
5. **Root cause** (from the source, `src/distribution/pareto.rs`): `cdf`
   computes `1.0 - (self.scale / x).powf(self.shape)` and `sf`
   `(self.scale / x).powf(self.shape)`. Cause 1 is the cancellation in
   `1.0 - ...`; cause 2 is the quotient `scale / x` leaving the normal range
   before the power brings it back. Confirmed with a copy that computes
   log(x_m / x) as -ln_1p((x - x_m) / x_m) for x <= 2 x_m (x - x_m is exact
   there) and `cdf` as `-(shape * log(x_m / x)).exp_m1()`; and `sf` as
   `(scale / x).powf(shape)` where the quotient is normal, otherwise
   `q.powf(shape) * (2^-1000).powf(shape)` with q = (scale 2^1000) / x (exact
   scaling), and exp(shape (ln scale - ln x)) only where q is still below
   2.2e-308. This copy passes every grid point (worst 2.2e-16 for `cdf`,
   5.6e-15 for `sf`) and all 20,000 random points (worst 4.4e-16 and
   7.9e-15).

## Notes

The grid's scales are powers of two, so x / x_m is exact and the cdf near
the scale can be computed to full precision from it. For other scales the
quotient is rounded, and the cdf's condition number in it is about
1 / (x / x_m - 1), so a careful method subtracts first: x - x_m is exact for
x <= 2 x_m. The random scan uses such scales.

Users cannot work around cause 1 with `sf`: near the scale, 1 - sf(x) is
the same subtraction. `pdf` and `ln_pdf` are not checked.
