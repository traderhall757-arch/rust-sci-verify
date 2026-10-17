# 042: statrs Gumbel cdf and sf return exactly 1 and 0 where x - location overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Gumbel`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location -2^1023 and
  scale 2^1023, `cdf(2^1023)` is 1 instead of 0.873 and `sf(2^1023)` is 0
  instead of 0.127), in a narrow range: location and x of opposite sign,
  both near the largest double. Off the grid also a small accuracy loss
  (up to 3.2e-13) in the far lower tail of `cdf`
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 4 of 1092 points of `reference/data/gumbel_cdf.jsonl`
  and 4 of 1092 of `reference/data/gumbel_sf.jsonl` (tolerance 1e-13). The
  exact lists are `GUMBEL_CDF` and `GUMBEL_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes:
  1. **x - mu above the largest double** (grid: mu = -2^1023, beta =
     2^1023, z = (x - mu) / beta = 2, 2.25, 2.5 and 2.75): the difference
     overflows to inf although z is moderate, so `cdf` returns exactly 1
     and `sf` exactly 0. Every other grid point passes (worst 4.8e-14 for
     `cdf`, 2.2e-16 for `sf`).
  2. **Rounding z in the steep lower tail** (off the grid only; the grid
     makes z exact on purpose, see Notes): F = exp(-exp(-z)) has condition
     number |z| exp(-z) in z, up to about 4,600 where F is still a normal
     number, so rounding (x - mu) / beta once costs up to about 5e-13.
     Seen at 102 of 20,000 random points, all with z from -6.55 to -5.22
     (F from 1.5e-305 to 9.7e-81), worst 3.1e-13; for example
     `Gumbel::new(0.0, 0.001)` gives `cdf(-0.0065)` 3.2e-13 off.

`Gumbel::inverse_cdf` fails elsewhere: [finding 043](043-statrs-gumbel-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Gumbel};

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let g = Gumbel::new(-big, big).unwrap(); // (location, scale)
    println!("{:?}", g.cdf(big)); // should be 0.8734230184931167 (z = 2)
    println!("{:?}", g.sf(big)); // should be 0.12657698150688335
    println!("{:?}", g.inverse_cdf(0.9)); // should be 1.1238883801828629e308
    let s = Gumbel::new(0.0, 1.0).unwrap();
    println!("{:?}", s.inverse_cdf(0.36787944117144233)); // should be 3.3784855259134224e-17
    println!("{:?}", s.inverse_cdf(0.36787944127144234)); // should be 2.7182823912192446e-10
    let t = Gumbel::new(0.0, 10.0).unwrap();
    println!("{:?}", t.inverse_cdf(0.3679)); // should be 0.0005588469011184362
    let n = Gumbel::new(0.0, 0.001).unwrap();
    println!("{:?}", n.cdf(-0.0065)); // should be 1.3572476073246106e-289 (off the grid)
    // Remedies: z = x / beta - mu / beta where x - mu overflows; log(-log p) as
    // log1p(-log1p(d)) with d = e p - 1 in double-double; beta (mu / beta - L).
    println!("{:?}", (-(-(big / big - -big / big)).exp()).exp()); // 0.8734230184931167
    let (p, e_lo) = (0.36787944117144233f64, 1.4456468917292502e-16);
    let prod = p * std::f64::consts::E;
    let d = (prod - 1.0) + (p.mul_add(std::f64::consts::E, -prod) + p * e_lo);
    println!("{:?}", -(-d.ln_1p()).ln_1p()); // 3.3784855259134224e-17
    let l = (-(0.9f64.ln())).ln();
    println!("{:?}", big * (-big / big - l)); // 1.1238883801828627e308: one ulp off
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below and the rows of finding 043.

## Expected

| Function (mu, beta, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| cdf (-2^1023, 2^1023, 2^1023) | 1.0 | 0.8734230184931167 | 1.0 | 0.14 |
| sf (-2^1023, 2^1023, 2^1023) | 0.0 | 0.12657698150688335 | 0.0 | 1 |
| cdf (-2^1023, 2^1023, 1.5 2^1023) | 1.0 | 0.9211936551755158 | 1.0 | 0.086 |
| sf (-2^1023, 2^1023, 1.5 2^1023) | 0.0 | 0.07880634482448418 | 0.0 | 1 |
| cdf (0, 0.001, -0.0065), off the grid | 1.3572476073250428e-289 | 1.3572476073246106e-289 | 1.3572476073250428e-289 | 3.2e-13 |

Second opinion: SciPy (`scipy.stats.gumbel_r`) computes (x - loc) / scale
the same way and is wrong at every row, so every row was also computed by
a second mpmath route: exp(-exp(-z)) and 1 - exp(-exp(-z)) at 1200 bits
(no expm1), which agree with the table to every printed digit. In
`reference/tests/test_functions.py` the generator's route is tested
against quadrature of the density exp(-t - exp(-t)) in t (both tails,
including the rows above), the series 1 - exp(-w) = sum of (-1)^(k+1)
w^k / k! in the upper tail (no expm1), and, for the standardisation only,
the exponential distribution of exp(-Z).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route (`_mp_gumbel_cdf` and
   `_mp_gumbel_sf`: exp(-w) and -expm1(-w) with w = exp(-(x - mu) / beta))
   is tested as described under "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, beta, x) (Python
   `random.Random(20261020)`), each run through `cdf` and `sf`: beta =
   10^v with v uniform in [-300, 300] and mu = 0 (a quarter) or +-10^u
   with u uniform in [-5, 5], both rounded to 6 significant digits (5% of
   the points instead have |mu| = 10^u, u in [300, 308.2], and beta = 10^v,
   v in [307, 308.2]); x = mu + beta z rounded to a double, with z uniform
   in [-7, 40] for half, 10^w (w uniform in [1.5, 3]) for a quarter and
   -10^w (w uniform in [-10, 0.85]) for a quarter; mpmath at 50 and 100
   digits (they agreed at all points). `cdf` fails 102 (cause 2, see
   above); `sf` fails none (worst 7.7e-14). A second scan
   (`random.Random(20261021)`) aims at cause 1: 2,000 points with
   mu = -10^u, beta = 10^v and x = 10^w, u, v and w uniform in
   [307.5, 308.25], 6 significant digits. `cdf` and `sf` each fail 831,
   every one where x - mu overflows (z from 1.08 to 9.57), returning
   exactly 1 and 0.
5. **Root cause** (from the source, `src/distribution/gumbel.rs`): `cdf`
   computes `(-(-(x - self.location) / self.scale).exp()).exp()` and `sf`
   `-(-(-(x - self.location) / self.scale).exp()).exp_m1()`. Cause 1 is
   `x - self.location` overflowing before the division brings it back;
   cause 2 is that the subtraction and the division each round z.
   Confirmed with a copy that computes z = x / beta - mu / beta where
   x - mu overflows, and otherwise carries z as z_hi + z_lo (an exact
   TwoSum for x - mu and a fused multiply-add for the division's
   remainder), then w_hi = exp(-z_hi), w_lo = -w_hi z_lo, F = exp(-w_hi)
   (1 - w_lo) and Q = -expm1(-(w_hi + w_lo)). This copy passes every grid point (worst
   4.8e-14 for `cdf`, 2.2e-16 for `sf`), all 20,000 random points (worst
   5.5e-14 and 2.4e-16) and all 2,000 targeted points (worst 2.1e-16 and
   1.2e-15).

## Notes

The grid's scales are powers of two and every x is mu + beta z exactly,
so z carries no rounding and the tolerance judges the evaluation alone
(the tolerance reason relies on this). Off the grid, cause 2 is an
accuracy loss rather than a wrong answer: a method that rounds z once
loses up to |z| exp(-z) 2^-53, about 5e-13, in the far lower tail, and only
a method that keeps z's rounding error (as the copy above does) stays
within the tolerance there. It is listed so that a fix covers it, and it
does not by itself raise the severity.
