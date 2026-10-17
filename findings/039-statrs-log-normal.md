# 039: statrs LogNormal cdf and sf have about 10 correct digits

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::LogNormal`,
  methods `cdf` and `sf` of `ContinuousCDF`)
- **Severity:** low (inaccuracy, up to 1.03e-10 relative, about 10 correct digits
  instead of 15 to 16)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 308 of 652 points of `reference/data/lognorm_cdf.jsonl`
  and 308 of 652 points of `reference/data/lognorm_sf.jsonl` (tolerance
  8e-13); the exact lists are `LOGNORMAL_CDF` and `LOGNORMAL_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  With w = (ln x - mu) / sigma, `cdf` fails at every grid point with
  -37.5 <= w <= -1 (4.9e-11 to 8.6e-11) and 1 <= w <= 2 (the same absolute
  error, diluted by a value near 1: 1.0e-12 at w = 2), for all four
  (mu, sigma) pairs; `sf` is the mirror image. |w| <= 0.5 passes, and so do
  the subnormal and zero results beyond w = -37.5 (relative to 2.2e-308).
  The quantile `inverse_cdf` is **not** affected: it passes all 548 untagged
  points of `lognorm_ppf` (worst 3.3e-14 as the report rounds it, tolerance 1e-13).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, LogNormal};

fn main() {
    let d = LogNormal::new(0.0, 1.0).unwrap();
    println!("{:?}", d.cdf(0.36787944117144233)); // 0.15865525394505725, should be 0.15865525393145705
    println!("{:?}", d.sf(1.9321599304402836e16)); // 4.605353009356527e-308, should be 4.605353009581957e-308
}
```

The "Crate" column below comes from a scratch program (not committed) that
reads (function, mu, sigma, x) lines from standard input and prints the
result with `{:?}`.

## Expected

| mu | sigma | x (w) | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|---|
| 0 | 1 | 0.36787944117144233 (-1) | `cdf` | 0.15865525394505725 | 0.15865525393145705 | 0.15865525393145707 | 8.6e-11 |
| 0 | 1 | 4.5399929762484854e-05 (-10) | `cdf` | 7.619853024559641e-24 | 7.61985302416053e-24 | 7.61985302416047e-24 | 5.2e-11 |
| 2 | 4 | 1.3336148155022614e-34 (-20) | `cdf` | 2.753624118411543e-89 | 2.753624118606234e-89 | 2.7536241186061556e-89 | 7.1e-11 |
| -3 | 0.5 | 3.581747930283181e-10 (-37.5) | `cdf` | 4.605353009356527e-308 | 4.605353009581979e-308 | 4.605353009581954e-308 | 4.9e-11 |
| -3 | 0.5 | 0.1353352832366127 (2) | `cdf` | 0.9772498680528374 | 0.9772498680518208 | 0.9772498680518208 | 1.0e-12 |
| 0 | 0.25 | 1.2840254166877414 (1) | `sf` | 0.15865525394505733 | 0.15865525393145713 | 0.15865525393145713 | 8.6e-11 |
| 0 | 1 | 1.9321599304402836e16 (37.5) | `sf` | 4.605353009356527e-308 | 4.605353009581957e-308 | 4.605353009581954e-308 | 4.9e-11 |

SciPy (`scipy.stats.lognorm` with `s = sigma`, `scale = exp(mu)`) agrees
with mpmath to 2.9e-14 or better at every row (second opinion).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode; it printed the two wrong values
   shown).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route Phi((ln x - mu) / sigma) through
   `_mp_erfc` is tested against direct quadrature of the density for the
   smaller tail (no erfc involved; 9 points out to |w| = 37.5, both tails,
   agreement to 47 digits or more) and cdf + sf = 1, in
   `reference/tests/test_functions.py`. A second test compares every grid
   point with `mpmath.ncdf`; that only checks the wiring (ln x, mu, sigma,
   the sign of w), because `ncdf` uses the same erfc formula for w < 0.
3. **Second opinion**: SciPy within tolerance at every row of "Expected".
4. **Off-grid**: 20,000 points with mu uniform in [-5, 5], sigma log-uniform
   in [0.1, 5] and w uniform in [-38.5, 38.5], x = exp(mu + sigma w) (Python
   `random.Random(20261018)`), each compared with mpmath at 50 digits. `cdf`
   fails 10,070 of them and `sf` 9,902 (tolerance 8e-13), worst 1.03e-10;
   failures span w = -37.62 to 2.10 for `cdf` and -2.10 to 37.63 for `sf`,
   matching the grid. 4,996 random quantiles (5,000 draws of p = U^1,
   U^10 or U^100 with U uniform, same seed, mu in [-5, 5], sigma in
   [0.1, 4]; the 4 draws that underflowed to 0 dropped) all pass, worst
   2.4e-14. The seed is only a number, the next in this project's series of
   date-shaped seeds.
5. **Root cause** (confirmed): `LogNormal::cdf` is `0.5 * erfc((location -
   ln x) / (scale * sqrt 2))` and `sf` is `0.5 * erfc((ln x - location) /
   (scale * sqrt 2))` (`src/distribution/log_normal.rs`), so they inherit the
   error of statrs's `erfc`, [finding 004](004-statrs-erfc.md), for
   |w / sqrt 2| >= 0.5 (the intervals of `erf_impl` with the mistyped
   constants), exactly as `Normal` does ([finding 009](009-statrs-normal.md)).
   A copy of statrs 0.19.1 with the 13 interval constants of `erf_impl`
   replaced by `<constant>_f32 as f64` (the fix of finding 004) removes every
   failure: 0 of the 40,000 off-grid checks and 0 of the 652 + 652 grid points
   fail, worst 4.5e-13 off-grid and 1.2e-13 on the grid. Fixing `erfc` fixes
   this finding.

## Notes

A finding of its own (rather than a citation of 004 or 009) for the reason
given in finding 009: users call `LogNormal::cdf`, not `erfc`, and the
affected range and size of the error are stated in its own variable w.
`pdf` and `ln_pdf` do not call `erfc`; not checked.
