# 041: statrs Pareto inverse_cdf returns inf where (1 - p)^(-1/shape) overflows but the quantile does not

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Pareto`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high (a silently wrong answer: with shape 2^-5 and scale
  2^-500, `inverse_cdf(0.9999999999)` is inf instead of 3.05e169), in a
  narrow range: shapes below about 0.05 with small scales and p close to 1
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 1 of 1215 points of `reference/data/pareto_ppf.jsonl`
  (tolerance 1e-13): shape 2^-5, scale 2^-500, p = 1 - 1e-10. The exact
  list is `PARETO_PPF` in `crates/check-statrs/src/expected.rs`, the report
  is [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). The
  quantile is x_m (1 - p)^(-1/alpha); the power alone exceeds 1.8e308
  wherever -log(1 - p) / alpha > 709.78, which needs alpha < 0.052 for
  p <= 1 - 2^-53. Where x_m is small enough the product is finite, but statrs
  returns inf. Every other grid point passes (worst 1.9e-14).

`Pareto::cdf` and `sf` fail for other reasons: [finding 040](040-statrs-pareto-cdf.md).

## Reproduction

This program was run in its own fresh project (the program of
[finding 040](040-statrs-pareto-cdf.md) prints the same line):

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, Pareto};

fn main() {
    let tiny = 2f64.powi(-500); // 3.054936363499605e-151
    let p = Pareto::new(tiny, 2f64.powi(-5)).unwrap(); // (scale, shape)
    println!("{:?}", p.inverse_cdf(0.9999999999)); // should be 3.054928274980468e169
    let h = (1.0f64 - 0.9999999999).powf(-0.5 * 32.0); // half the power, applied twice
    println!("{:?}", tiny * h * h); // 3.054928274980469e169: one ulp from the true value
}
```

Run in release and debug mode (same output): `inf`, then the remedy's value.

## Expected

| alpha, x_m, p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 2^-5, 2^-500, 0.9999999999 | inf | 3.054928274980468e169 | inf | inf |

Second opinion: SciPy (`scipy.stats.pareto.ppf`) also returns inf, so the
value was also computed as x_m (1 - p)^-32 directly at 1200 bits
(3.0549282749804682e169). In `reference/tests/test_functions.py` the
generator's route, x_m exp(-log1p(-p) / alpha), is tested against that
power and against cdf(ppf(p)) = p.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: mpmath at 50 and 100 digits round to the same double
   (the generator enforces this); the route is tested as above.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (alpha, x_m, p) (Python
   `random.Random(20261019)`, after the 20,000 points of finding 040):
   alpha = 10^u, u uniform in [-1.6, 2], and x_m = 10^v, v uniform in
   [-300, 300], both rounded to 6 significant digits; p = 10^w with w
   uniform in [-320, 0] for half, uniform in (0, 1) for a quarter and
   1 - 10^w with w uniform in [-16, 0] for a quarter (a p that rounds to 0
   or 1 is replaced by 1/2); mpmath at 50 and 100
   digits (they agreed at all points). `inverse_cdf` fails 44, every one
   returning inf, with alpha <= 0.0474 and p >= 0.99999999258. No other
   failure: the worst error at the other points is 6.1e-14.
5. **Root cause** (from the source, `src/distribution/pareto.rs`):
   `inverse_cdf` computes `self.scale * (1.0 - p).powf(-1.0 / self.shape)`;
   the power overflows before the multiplication by the scale. Confirmed
   with a copy that, where that power is infinite and p < 1, computes
   h = (1 - p)^(-1/(2 alpha)) and returns `scale * h * h`: it passes every
   grid point (worst 1.9e-14) and all 20,000 random points (worst 9.7e-14).

## Notes

1 - p is exact for p >= 1/2, so the overflowing case loses nothing in the
subtraction. For p below 1/2 statrs's 1 - p is rounded, which costs up to
2^-53 / alpha in the result (3.6e-15 at shape 2^-5): within the tolerance.

The tolerance reason's bound (345.4 2^-53 = 3.8e-14 for rounding 1/alpha)
holds for this grid. Off the grid, |log(x / x_m)| can reach about 1,400
where 1/alpha is not a double, and rounding 1/alpha then costs up to about
1.6e-13 for any method that uses it: a future off-grid failure of that size
would not by itself be a bug.
