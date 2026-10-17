# 047: statrs Cauchy inverse_cdf is finite at p = 0 and 1, loses accuracy next to them, and returns inf where scale times the tangent overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Cauchy`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location 0 and scale 1,
  `inverse_cdf(0.0)` is -1.6e16 instead of -inf, `inverse_cdf(1e-20)` is
  -1.6e16 instead of -3.2e19 and `inverse_cdf(1e-10)` is off by 5e-7)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 62 of 186 points of `reference/data/cauchy_ppf.jsonl`
  (tolerance 1e-13). The exact list is `CAUCHY_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). `inverse_cdf`
  computes `self.location + self.scale * (f64_consts::PI * (x - 0.5)).tan()`.
  Four causes:
  1. **Finite at p = 0 and p = 1** (10 grid points, tagged pole): pi / 2
     is not a double, so tan(pi (0 - 0.5)) is -1.633123935319537e16, not
     -inf, and the result is mu - 1.6e16 b (and mu + 1.6e16 b at p = 1).
     Only with scale 2^1023 does the product overflow to the right
     infinity. Rule 4: a finite value at a pole is crate wrong.
  2. **Rounding next to p = 0 and p = 1** (50 grid points: p from 5e-324
     to 1e-5 and from 1 - 1e-5 to 1 - 1e-15, at every location and scale;
     5 of them are tagged overflow because the true quantile is beyond the
     largest double, and statrs returns a finite number there). The
     quantile is about -b / (pi p), set by the distance pi p between
     pi (p - 1/2) and -pi / 2; but `x - 0.5` and the product with pi each
     round with an absolute error of about 1e-16, so the relative error is
     about 1e-16 / p: 5e-7 at p = 1e-10. For p below 2^-54, `x - 0.5`
     rounds to exactly -0.5 and the result is the p = 0 value, mu - 1.6e16 b,
     whatever p is (-1.6e16 instead of -3.2e19 at p = 1e-20). Near p = 1,
     1 - p is lost the same way.
  3. **Not 0 where the quantile is 0** (1 grid point: location -2^1023,
     scale 2^1023, p = 3/4): tan of the double nearest pi / 4 is
     0.9999999999999999, not 1, so mu + b T is -9.98e291 instead of 0.
  4. **b T above the largest double** (1 grid point: location -2^1023,
     scale 2^1023, p = 7/8): `self.scale * tan(..)` overflows to inf
     although mu + b T = 1.27e308 is finite (cause 2 of finding 045 in
     `Laplace`).
  Every other point passes (worst 3.1e-14, at location -10, scale 1/4,
  p = 0.999, where cause 2 has begun); the 22 points where the true
  quantile is beyond the largest double and statrs returns the right
  infinity also pass.

`Cauchy::cdf` and `sf` fail for other reasons: [finding 046](046-statrs-cauchy-cdf.md).
Cause 4 is the same as cause 2 of [finding 045](045-statrs-laplace-inverse-cdf.md)
(`Laplace`) and of [finding 043](043-statrs-gumbel-inverse-cdf.md) (`Gumbel`).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Cauchy, ContinuousCDF};
use std::f64::consts::{FRAC_1_PI, PI};

/// z = (x - mu) / b, with x / b - mu / b where x - mu overflows.
fn z_of(mu: f64, b: f64, x: f64) -> f64 {
    let d = x - mu;
    if d.is_finite() { d / b } else { x / b - mu / b }
}

/// F(z) without cancellation: atan(1 / |z|) / pi in the lower tail.
fn cdf_z(z: f64) -> f64 {
    if z < -1.0 {
        (-1.0 / z).atan() / PI
    } else if z > 1.0 {
        1.0 - (1.0 / z).atan() / PI
    } else {
        0.5 + z.atan() / PI
    }
}

/// b / tan(pi q) for 0 <= q < 1/4 (q exact), as b / (pi q) for tiny q.
fn b_cot(b: f64, q: f64) -> f64 {
    if q < 2f64.powi(-30) {
        let s = b * FRAC_1_PI;
        if s >= f64::MIN_POSITIVE { s / q } else { b / q * FRAC_1_PI }
    } else {
        b / (PI * q).tan()
    }
}

/// tan(pi d) for |d| <= 1/4 (d exact), exactly +-1 at d = +-1/4.
fn tan_pi(d: f64) -> f64 {
    let a = d.abs();
    let t = if a <= 0.125 {
        (PI * a).tan()
    } else {
        let u = (PI * (0.25 - a)).tan(); // 1/4 - a is exact
        (1.0 - u) / (1.0 + u)
    };
    t.copysign(d)
}

fn ppf_fix(mu: f64, b: f64, p: f64) -> f64 {
    let bt = if p < 0.25 {
        -b_cot(b, p)
    } else if p > 0.75 {
        b_cot(b, 1.0 - p)
    } else {
        b * tan_pi(p - 0.5)
    };
    let x = mu + bt;
    if bt.is_infinite() && x.is_infinite() {
        // b T overflowed; mu + b T may still be finite.
        let t = if p < 0.25 { -b_cot(1.0, p) } else if p > 0.75 { b_cot(1.0, 1.0 - p) } else { tan_pi(p - 0.5) };
        if t.is_finite() { b * (mu / b + t) } else { x }
    } else {
        x
    }
}

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let s = Cauchy::new(0.0, 1.0).unwrap(); // (location, scale)
    println!("{:?}", s.cdf(-1e16)); // should be 3.183098861837907e-17
    println!("{:?}", s.cdf(-1e10)); // should be 3.1830988618379065e-11
    println!("{:?}", s.sf(1e10)); // should be 3.1830988618379065e-11
    let l = Cauchy::new(-big, big).unwrap();
    println!("{:?}", l.cdf(big)); // should be 0.8524163823495667 (z = 2)
    println!("{:?}", l.sf(big)); // should be 0.14758361765043326
    println!("{:?}", s.inverse_cdf(0.0)); // should be -inf
    println!("{:?}", s.inverse_cdf(1.0)); // should be inf
    println!("{:?}", s.inverse_cdf(1e-20)); // should be -3.1830988618379067e19
    println!("{:?}", s.inverse_cdf(1e-10)); // should be -3183098861.837907
    println!("{:?}", l.inverse_cdf(0.75)); // should be 0
    println!("{:?}", l.inverse_cdf(0.875)); // should be 1.2711610061536464e308
    // Remedies: atan(1 / |z|) / pi in the tails, x / b - mu / b where x - mu
    // overflows; -b / tan(pi p) (b / (pi p) for tiny p), tan(pi / 4 - y) next to
    // p = 3/4, and b (mu / b + T) where b T overflows.
    println!("{:?}", cdf_z(z_of(0.0, 1.0, -1e16)));
    println!("{:?}", cdf_z(z_of(0.0, 1.0, -1e10)));
    println!("{:?}", cdf_z(-z_of(0.0, 1.0, 1e10)));
    println!("{:?}", cdf_z(z_of(-big, big, big)));
    println!("{:?}", cdf_z(-z_of(-big, big, big)));
    println!("{:?}", ppf_fix(0.0, 1.0, 0.0));
    println!("{:?}", ppf_fix(0.0, 1.0, 1.0));
    println!("{:?}", ppf_fix(0.0, 1.0, 1e-20));
    println!("{:?}", ppf_fix(0.0, 1.0, 1e-10));
    println!("{:?}", ppf_fix(-big, big, 0.75));
    println!("{:?}", ppf_fix(-big, big, 0.875));
}
```

This is the program that was run (release and debug mode give the same
output). It prints the rows of finding 046 (0.0, 3.183098229442294e-11,
3.183098229442294e-11, 1.0, 0.0), then -1.633123935319537e16,
1.633123935319537e16, -1.633123935319537e16, -3183097229.936138,
-9.9792015476736e291 and inf for the rows below, then the remedies'
values: 3.183098861837907e-17, 3.183098861837907e-11 (twice),
0.8524163823495667, 0.14758361765043326, -inf, inf,
-3.183098861837907e19, -3183098861.837907, 0.0 and
1.2711610061536462e308, each within an ulp of the true value.

## Expected

| (mu, b, p) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (0, 1, 0) | -1.633123935319537e16 | -inf (pole) | -inf | pole |
| (0, 1, 1) | 1.633123935319537e16 | +inf (pole) | inf | pole |
| (0, 1, 1e-20) | -1.633123935319537e16 | -3.1830988618379067e19 | -3.183098861837907e19 | 0.9995 |
| (0, 1, 1e-10) | -3183097229.936138 | -3183098861.837907 | -3183098861.8379064 | 5.1e-7 |
| (0, 1, 1 - 1e-15) | 315471663322406.1 | 318564507734592.1 | 318564507734592.2 | 9.7e-3 |
| (-2^1023, 2^1023, 3/4) | -9.9792015476736e291 | 0 | 1.99584030953472e292 | inf |
| (-2^1023, 2^1023, 7/8) | inf | 1.2711610061536464e308 | inf | inf |

Second opinion: SciPy (`scipy.stats.cauchy`) agrees with mpmath at the
cause-1 and cause-2 rows above (it is wrong at p = 5e-324 with scale
2^-500, a cause-2 point, where statrs is also wrong), and is wrong at the
cause-3 and cause-4 rows, so those were checked by hand: at p = 3/4,
tan(pi / 4) = 1 and mu + b = 0 exactly; at p = 7/8, tan(3 pi / 8) =
1 + sqrt 2, so the quantile is sqrt 2 times 2^1023 = 1.2711610061536464e308.
The generator's route (mu - b cospi(p) / sinpi(p), no subtraction from p)
is tested in `reference/tests/test_functions.py` against the root of
cdf(x) = p at 1200 bits (including p = 5e-324, 1 - 1e-15 and 7/8, and the
scale-2^1023 parameters).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: mpmath at 50 and 100 digits round to the same double
   (the generator enforces this); the route is tested as above. p = 0 and 1
   are tagged pole by the rule in `_cauchy_ppf_tags`, overflow by the
   generator's rule.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, b, p) (Python `random.Random(20261025)`,
   after the 20,000 points of finding 046), with mu and b drawn as there;
   p = 10^w with w uniform in [-320, 0] for a quarter, uniform in (0, 1)
   for a quarter, 1 - 10^w with w uniform in [-16, 0] for a quarter and
   1/2 +- 10^w with w uniform in [-17, -1] for a quarter (a p that rounds
   to 0 or 1 is replaced by 1/2); mpmath at 50 and 100 digits (they agreed
   at all points). `inverse_cdf` fails 6,075: 6,069 with min(p, 1 - p) at
   most 4.7e-4 (cause 2; relative error from 1.0e-13 at p near 1e-3 to 1)
   and 6 returning inf where the quantile is finite (cause 4). The two
   overflow scans of finding 046 (`random.Random(20261026)` and
   `random.Random(20261027)`, interleaved with its points) also draw
   p uniform in (0, 1) for each of their 2,000 (mu, b): 67 and 68
   failures, every one returning +-inf where the quantile is finite
   (cause 4). Causes 1 and 3 are at fixed inputs and need no scan.
5. **Root cause** (from the source, `src/distribution/cauchy.rs`):
   `inverse_cdf` checks `(0.0..=1.0).contains(&x)` and then computes
   `self.location + self.scale * (f64_consts::PI * (x - 0.5)).tan()`.
   Causes 1 to 3 come from evaluating the tangent at a rounded multiple of
   the rounded pi, cause 4 from the product overflowing. Confirmed with
   the copy `ppf_fix` above (cotangent of pi p or pi (1 - p) away from
   the middle, b / (pi p) for tiny p, tan(pi / 4 - y) next to p = 1/4 and
   3/4, b (mu / b + T) where b T overflows): it passes every grid point
   and all 20,000 + 2,000 + 2,000 random points (worst 6.0e-15 on the
   first scan, 3.3e-14 on the second and third, from mu + b T
   cancelling).
