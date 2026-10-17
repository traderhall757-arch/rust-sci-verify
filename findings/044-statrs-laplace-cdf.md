# 044: statrs Laplace cdf and sf return 1 and 0 where x - location overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Laplace`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location -2^1023 and
  scale 2^1023, `cdf(2^1023)` is 1 instead of 0.932 and `sf(2^1023)` is 0
  instead of 0.068)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 4 of 1828 points of `reference/data/laplace_cdf.jsonl`
  and the same 4 of `reference/data/laplace_sf.jsonl` (tolerance 1e-13).
  The exact lists are `LAPLACE_CDF` and `LAPLACE_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). All 4 have
  location -2^1023, scale 2^1023 and x = 2^1023 to 1.75 2^1023 (z = (x -
  location) / scale from 2 to 2.75): `x - self.location` overflows to inf
  although z is moderate, so `cdf` returns exactly 1 and `sf` exactly 0.
  Off the grid this happens wherever |x - location| exceeds the largest
  double (about 1.8e308), for any scale (see the checklist); where
  x - location is below -1.8e308 the results are 0 and 1 instead. Every other
  grid point passes (worst 2.2e-16 for both).

The quantile fails for other reasons: [finding 045](045-statrs-laplace-inverse-cdf.md).
The same overflow in `Gumbel` is [finding 042](042-statrs-gumbel-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Laplace};

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let l = Laplace::new(-big, big).unwrap(); // (location, scale)
    println!("{:?}", l.cdf(big)); // should be 0.9323323583816937 (z = 2)
    println!("{:?}", l.sf(big)); // should be 0.06766764161830635
    println!("{:?}", l.cdf(1.5 * big)); // should be 0.9589575006880506 (z = 2.5)
    println!("{:?}", l.sf(1.5 * big)); // should be 0.0410424993119494
    println!("{:?}", l.inverse_cdf(0.95)); // should be 1.1708241396246928e308
    let m = Laplace::new(4.67591e307, 1.50284e308).unwrap();
    println!("{:?}", m.inverse_cdf(0.14745759143266557)); // should be -1.367478051182768e308
    let s = Laplace::new(0.0, 1.0).unwrap();
    for p in [0.0, 1.0] {
        // should be -inf and inf
        println!("{:?}", std::panic::catch_unwind(|| s.inverse_cdf(p)));
    }
    // Remedies: z = x / b - mu / b where x - mu overflows; b (mu / b + L) where
    // b L overflows.
    println!("{:?}", 1.0 - (-(big / big - -big / big)).exp() / 2.0); // 0.9323323583816936
    println!("{:?}", big * (-big / big - (2.0 - 2.0 * 0.95f64).ln())); // 1.170824139624693e308
}
```

This is the program that was run (release and debug mode give the same
output). It prints 1.0, 0.0, 1.0, 0.0 for the four rows below, then the
rows of finding 045 (inf, -inf and two panics with the message
"p must be in [0, 1]"), then the remedies' values, each one ulp from the
true value.

## Expected

| Function (mu, b, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| cdf (-2^1023, 2^1023, 2^1023) | 1.0 | 0.9323323583816937 | 1.0 | 0.073 |
| sf (-2^1023, 2^1023, 2^1023) | 0.0 | 0.06766764161830635 | 0.0 | 1 |
| cdf (-2^1023, 2^1023, 1.5 2^1023) | 1.0 | 0.9589575006880506 | 1.0 | 0.043 |
| sf (-2^1023, 2^1023, 1.5 2^1023) | 0.0 | 0.0410424993119494 | 0.0 | 1 |

Second opinion: SciPy (`scipy.stats.laplace`) computes (x - loc) / scale
the same way and is wrong at every row, so every row was also computed by
a second mpmath route: e^-|z| / 2 and 1 - e^-|z| / 2 at 1200 bits, which
agree with the table to every printed digit. In
`reference/tests/test_functions.py` the generator's route is tested
against quadrature of the density e^-|t| / 2 in t (both tails, including
the first and last rows' parameters), and against half the exponential
survival function for the branches and the mirror Q(x) = F(2 mu - x).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route (`_mp_laplace_cdf`: e^z / 2 for
   z < 0 and 1 - e^-z / 2 for z >= 0; `_mp_laplace_sf` by symmetry) is
   tested as described under "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, b, x) (Python `random.Random(20261022)`),
   each run through `cdf` and `sf`: b = 10^v with v uniform in [-300, 300]
   and mu = 0 (a quarter) or +-10^u with u uniform in [-5, 5], both rounded
   to 6 significant digits (5% of the points instead have |mu| = 10^u,
   u in [300, 308.2], and b = 10^v, v in [307, 308.2]); x = mu + b z
   rounded to a double, with z uniform in [-40, 40] for half and
   +-10^w for the rest (w uniform in [1.5, 2.9] for a quarter, in
   [-10, 1.5] for a quarter), and x = mu where mu + b z overflows; mpmath at 50 and 100 digits (they agreed at
   all points). `cdf` and `sf` fail none. A second scan
   (`random.Random(20261023)`) aims at the overflow: 2,000 points with
   mu = -10^u, b = 10^v and x = 10^w, u, v and w uniform in
   [307.5, 308.25], 6 significant digits. `cdf` and `sf` each fail 813,
   every one where x - mu overflows (|z| from 1.05 to 10.4), returning
   exactly 1 and 0. A third scan (`random.Random(20261024)`) covers the
   mirror case: 2,000 points with mu = 10^u, b = 10^v and x = -10^w, drawn
   as in the second. `cdf` and `sf` each fail 819, every one where x - mu
   overflows to -inf (|z| from 1.05 to 9.42), returning exactly 0 and 1.
5. **Root cause** (from the source, `src/distribution/laplace.rs`): `cdf`
   and `sf` compute `(-(x - self.location).abs() / self.scale).exp() / 2.`
   and then `1. - y` or `y`. `x - self.location` overflows to inf before
   the division brings it back, so y is 0. Confirmed with a copy that
   computes z = x / b - mu / b where x - mu overflows (and otherwise
   carries |z| as hi + lo, from an exact TwoSum for x - mu and a fused
   multiply-add for the division's remainder, with y = e^-hi (1 - lo) / 2).
   This copy passes every grid point and all 40,000 + 4,000 + 4,000
   random evaluations (worst 3.1e-16 on the first scan, 1.4e-15 on the
   second, 1.3e-15 on the third).

## Notes

Unlike `Gumbel` (finding 042's notes), rounding z costs little here: F's
condition number in z is |z| for z < 0, and the first scan reaches |z| of
about 800, so each rounding of z costs at most about 9e-14. statrs's own z
(one rounding in the subtraction, one in the division) stays within the
tolerance at every random point (worst 6.2e-14 for `cdf`, 5.5e-14 for
`sf`).
