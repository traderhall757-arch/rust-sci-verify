# 049: statrs Levy inverse_cdf returns inf where scale / (2 erfc_inv(p)^2) overflows although the quantile is finite

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Levy`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high, narrow range (silently wrong answers: with location
  -2^1023 and scale 2^1023, `inverse_cdf(0.5)` is inf instead of 1.077e308;
  only where the location is near minus the largest double)
- **Status:** approved (Claude, checklist complete, 2026-10-17)
- **Affected inputs:** 5 of 180 points of `reference/data/levy_ppf.jsonl`
  (tolerance 1e-13). The exact list is `LEVY_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). `inverse_cdf`
  computes `self.mu + 0.5 * self.c / (erfc_inv(x).powf(2.0))`. With
  location -2^1023 and scale 2^1023, 0.5 c / y^2 (y = erfc_inv(p))
  overflows to inf for every p above erfc(1/2) = 0.4795, while the
  quantile mu + 0.5 c / y^2 stays finite up to p = erfc(1/sqrt(6)) =
  0.5637. The grid points p = 1/2 - 1e-10, 1/2 - 2^-54, 1/2,
  1/2 + 2^-53 and 1/2 + 1e-10 fall in that interval; statrs returns inf
  at all five (true values 1.0769164652e308 to 1.0769164671e308).
  Every other point passes (worst 9.1e-16, at p = 1/4 with the same
  parameters), including p = 0 (exactly the location) and p = 1 (inf,
  tagged pole), and the 10 points tagged overflow, where statrs's inf is
  right.

The cdf and sf fail for other reasons: [finding 048](048-statrs-levy-cdf.md).
The same overflow of an intermediate in `Cauchy` and `Laplace` is
[finding 047](047-statrs-cauchy-inverse-cdf.md) and
[finding 045](045-statrs-laplace-inverse-cdf.md).

## Reproduction

The program of [finding 048](048-statrs-levy-cdf.md) (its last two lines):

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Levy};

fn main() {
    let big = 2f64.powi(1023); // 8.98846567431158e307
    let l = Levy::new(-big, big).unwrap(); // (location, scale)
    println!("{:?}", l.inverse_cdf(0.5)); // should be 1.0769164661541097e308
    println!("{:?}", l.inverse_cdf(0.4999999999)); // should be 1.0769164652323077e308
}
```

Run in release and debug mode (as part of finding 048's program) with the
same output: inf and inf.

## Expected

| Function (mu, c, p) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| inverse_cdf (-2^1023, 2^1023, 1/2) | inf | 1.0769164661541097e308 | inf | inf |
| inverse_cdf (-2^1023, 2^1023, 1/2 - 1e-10) | inf | 1.0769164652323077e308 | inf | inf |
| inverse_cdf (-2^1023, 2^1023, 1/4) | -2.1960221141866992e307 | -2.1960221141867012e307 | -2.1960221141866982e307 | 9.1e-16 (passes) |

Second opinion: SciPy (`scipy.stats.levy.ppf`) also returns inf at every
failing row (it forms scale / z^2 the same way). The generator's route,
mu + c / z^2 with z = Phi^-1(p / 2) at 60 extra bits, is tested in
`reference/tests/test_functions.py` against the root of the generator's cdf
(erfc, no inverse function) found by the secant method at 1200 bits, at
p = 1/2 - 2^-54 and 1/2 + 2^-53 with these parameters among others.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (finding 048's program, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the route is tested as described under
   "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: a scratch copy of statrs's `inverse_cdf` (with and without
   the remedy) and mpmath at 50 and 100 digits (they agreed at all points).
   Scan 4 (Python `random.Random(20261031)`): 5,000 random (mu, c, p):
   c = 10^v with v uniform in [-300, 300] and mu = 0 (a quarter) or
   +-10^u with u uniform in [-5, 5], 6 significant digits (10% of the
   points instead have mu = -10^u, c = 10^v, u and v uniform in
   [307, 308.2]); p = 10^-w with w uniform in [0, 300] for half,
   1 - 10^-w with w uniform in [0, 16] for a quarter, uniform in (0, 1)
   for a quarter. statrs fails 5, every one where 0.5 c / y^2 overflows
   while the quantile is finite. Scan 5 (`random.Random(20261032)`):
   2,000 points with mu = -10^u, c = 10^v (u, v uniform in [307, 308.25],
   6 significant digits) and p uniform in (0, 1): statrs fails 79; 78 are
   this overflow (p from 0.338 to 0.815, depending on mu and c). The 79th
   (mu = -8.536e307, c = 8.43532e307, p = 0.31991...) fails in the
   remedied copy too, by 1.06e-13: there mu + 0.5 c / y^2 cancels (the
   quantile is -9.4e304, 1815 times smaller than |mu|), which costs any
   method that rounds 0.5 c / y^2 up to about 1815 ulps. That is not a
   statrs error, and the tolerance does not claim otherwise: on the grid
   this ratio is at most 7.2 (tested).
5. **Root cause** (from the source, `src/distribution/levy.rs`): the
   intermediate `0.5 * self.c / (erfc_inv(x).powf(2.0))` overflows before
   the location is added. Confirmed with a copy that computes
   r = mu + c q with q = 0.5 / y^2 and, where r is not finite, (mu / c + q) c
   instead. This copy passes every grid point (worst 9.1e-16) and every
   point of both scans apart from the cancelling point above (worst
   passing 7.0e-15 on scan 4, 6.8e-14 on scan 5, both where mu and c q
   nearly cancel).

## Notes

statrs's `erfc_inv` itself is accurate on this grid: away from the five
points every quantile is within 9.1e-16, also next to p = 0 (p = 5e-324)
and p = 1 (1 - 2^-53).
