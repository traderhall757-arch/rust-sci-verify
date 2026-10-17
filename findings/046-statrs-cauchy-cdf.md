# 046: statrs Cauchy cdf and sf lose all digits in the tails, and return 1 and 0 where x - location overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Cauchy`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location 0 and scale 1,
  `cdf(-1e16)` is 0 instead of 3.2e-17 and `cdf(-1e10)` is off by 2e-7;
  with location -2^1023 and scale 2^1023, `cdf(2^1023)` is 1 instead of
  0.852)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 56 of 1762 points of `reference/data/cauchy_cdf.jsonl`
  and 56 of `reference/data/cauchy_sf.jsonl` (tolerance 1e-13). The exact
  lists are `CAUCHY_CDF` and `CAUCHY_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes:
  1. **Cancellation in the tail** (52 points each): `cdf` computes
     atan(z) / pi + 0.5 with z = (x - location) / scale. In the lower
     tail atan(z) / pi is close to -1/2, so the sum keeps only the last
     bits of a number near 1/2: the true value is about 1 / (pi |z|), and
     about log10 |z| digits are lost. The grid fails from |z| = 1000
     (relative error 1.3e-13) and every point from |z| = 1e16 on is
     exactly 0 (true values down to 3.2e-309). `sf` is the mirror image
     in the upper tail. Five of the six (location, scale) pairs are affected;
     the sixth (-2^1023, 2^1023) has no tail points on the grid (|z| <= 2.75).
  2. **x - location above the largest double** (4 points each: location
     -2^1023, scale 2^1023, x = 2^1023 to 1.75 2^1023, so z = 2 to 2.75):
     `x - self.location` overflows to inf although z is moderate, so `cdf`
     returns exactly 1 and `sf` exactly 0 (finding 044's cause in
     `Laplace`).
  Every other point passes (worst 2.1e-14 for both, at z = -100 and 100,
  where the cancellation has begun).

The quantile fails for other reasons: [finding 047](047-statrs-cauchy-inverse-cdf.md).
The same overflow in `Laplace` and `Gumbel` is [finding 044](044-statrs-laplace-cdf.md)
and [finding 042](042-statrs-gumbel-cdf.md); the same `1 - (something
near 1)` pattern is behind findings 023 and 040.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Cauchy, ContinuousCDF};

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let s = Cauchy::new(0.0, 1.0).unwrap(); // (location, scale)
    println!("{:?}", s.cdf(-1e16)); // should be 3.183098861837907e-17
    println!("{:?}", s.cdf(-1e10)); // should be 3.1830988618379065e-11
    println!("{:?}", s.sf(1e10)); // should be 3.1830988618379065e-11
    let l = Cauchy::new(-big, big).unwrap();
    println!("{:?}", l.cdf(big)); // should be 0.8524163823495667 (z = 2)
    println!("{:?}", l.sf(big)); // should be 0.14758361765043326
    // ... then the rows of finding 047 and the remedies (see there).
}
```

This program (the first part of the one in finding 047) was run in
release and debug mode with the same output: 0.0, 3.183098229442294e-11,
3.183098229442294e-11, 1.0 and 0.0.

## Expected

| Function (mu, b, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| cdf (0, 1, -1e3) | 0.0003183097800805168 | 0.0003183097800805589 | 0.0003183097800805589 | 1.3e-13 |
| cdf (0, 1, -1e10) | 3.183098229442294e-11 | 3.1830988618379065e-11 | 3.183098861837907e-11 | 2.0e-7 |
| cdf (0, 1, -1e15) | 2.7755575615628914e-16 | 3.1830988618379067e-16 | 3.183098861837907e-16 | 0.13 |
| cdf (0, 1, -1e16) | 0.0 | 3.183098861837907e-17 | 3.183098861837907e-17 | 1 |
| cdf (0, 1, -1e308) | 0.0 | 3.183098861837906e-309 | 3.183098861837906e-309 | 1 |
| cdf (-2^1023, 2^1023, 2^1023) | 1.0 | 0.8524163823495667 | 1.0 | 0.17 |
| sf (-2^1023, 2^1023, 2^1023) | 0.0 | 0.14758361765043326 | 0.0 | 1 |

`sf` gives the same numbers at the mirrored points (x replaced by
2 mu - x).

Second opinion: SciPy (`scipy.stats.cauchy`) agrees with mpmath within an
ulp in the tails (all cause-1 rows), but computes (x - loc) / scale the same way as
statrs and is wrong at every cause-2 row. The generator's route
(acot(-z) / pi for z < 0, 1 - acot(z) / pi otherwise, so no cancellation)
is tested in `reference/tests/test_functions.py` against quadrature of the
density 1 / (pi (1 + t^2)) (both tails, including -1e308, the scale-2^-500
tail and the cause-2 parameters) and against the series of acot for the
branches and the mirror Q(x) = F(2 mu - x).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program of finding 047, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the route is tested as described under
   "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, b, x) (Python `random.Random(20261025)`),
   each run through `cdf` and `sf`: b = 10^v with v uniform in [-300, 300]
   and mu = 0 (a quarter) or +-10^u with u uniform in [-5, 5], both rounded
   to 6 significant digits (5% of the points instead have |mu| = 10^u,
   u in [300, 308.2], and b = 10^v, v in [307, 308.2]); x = mu + b z
   rounded to a double, with z uniform in [-40, 40] for half,
   +-10^w with w uniform in [1.5, 300] for a quarter and in [-10, 1.5] for
   a quarter, and x = mu where mu + b z overflows; mpmath at 50 and 100
   digits (they agreed at all points). `cdf` fails 1,358 and `sf` 1,405,
   every one in its own tail (z below -1696 for `cdf`, above 2318 for
   `sf`; cause 1). A second scan (`random.Random(20261026)`) aims at the
   overflow: 2,000 points with mu = -10^u, b = 10^v and x = 10^w, u, v
   and w uniform in [307.5, 308.25], 6 significant digits. `cdf` and `sf`
   each fail 803, every one where x - mu overflows (|z| from 1.02 to
   9.86), returning exactly 1 and 0 (cause 2). A third scan
   (`random.Random(20261027)`) covers the mirror case (mu = 10^u,
   x = -10^w): 825 failures each, every one where x - mu overflows to
   -inf, returning exactly 0 and 1.
5. **Root cause** (from the source, `src/distribution/cauchy.rs`): `cdf` is
   `(1.0 / f64_consts::PI) * ((x - self.location) / self.scale).atan() +
   0.5` and `sf` the same with `self.location - x`. The final `+ 0.5`
   cancels (cause 1), and `x - self.location` overflows before the
   division brings it back (cause 2). Confirmed with a copy that computes
   z = x / b - mu / b where x - mu overflows, and F = atan(-1 / z) / pi
   for z < -1, 1 - atan(1 / z) / pi for z > 1 and 0.5 + atan(z) / pi
   otherwise (sf as F(-z)). This copy passes every grid point and all
   20,000 + 2,000 + 2,000 random points of each scan (worst 4.2e-16 on
   the first, 4.2e-16 on the second, 4.1e-16 on the third).

## Notes

Rounding z costs little: F's condition number in z, |z F'(z) / F(z)|, is
at most 1, so statrs's two roundings (subtraction, division) are harmless
by themselves; every cause-1 failure comes from the final addition.
