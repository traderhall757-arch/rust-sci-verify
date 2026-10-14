# 004: statrs erfc has about 10 correct digits

- **Crate / version:** statrs 0.19.1 (`statrs::function::erf::erfc`)
- **Severity:** low (inaccuracy, up to 1.0e-10 relative, about 10 correct digits
  instead of 15 to 16)
- **Status:** approved by owner (2026-10-14)
- **Affected inputs:** 111 of 281 points of `reference/data/erfc.jsonl`
  (tolerance 1e-13); the exact list is `ERFC` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  Every grid point with 0.5 <= x <= 26.5 fails (9e-12 to 1.0e-10), and so do
  x = -1.75 to -0.5, where erfc(x) = 2 - erfc(-x) carries the error at a smaller
  relative size (2.9e-13 to 3.2e-11). Beyond x = 26.5 the results are subnormal and
  within the error metric. SciPy's worst error on the same grid is 4.5e-16,
  apart from two subnormal results it flushes to zero.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::function::erf::erfc;

fn main() {
    println!("{:?}", erfc(0.5)); // 0.4795001222363462, should be 0.4795001221869535
    println!("{:?}", erfc(2.25)); // 0.0014627165868018256, should be 0.0014627165866811518
}
```

## Expected

| Input x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 0.5 | 0.4795001222363462 | 0.4795001221869535 | 0.4795001221869535 | 1.0e-10 |
| 2.25 | 0.0014627165868018256 | 0.0014627165866811518 | 0.0014627165866811515 | 8.2e-11 |
| 10 | 2.0884875837105948e-45 | 2.088487583762545e-45 | 2.0884875837625446e-45 | 2.5e-11 |
| 20 | 5.3958656116575216e-176 | 5.395865611607901e-176 | 5.395865611607901e-176 | 9.2e-12 |
| -1 | 1.842700792942715 | 1.8427007929497148 | 1.8427007929497148 | 3.8e-12 |

## Notes

The true values come from `mpmath.erfc` (and, for x >= 1e6, the asymptotic
series of erfcx), tested against erfc(x) = Gamma(1/2, x^2)/sqrt(pi) (DLMF 7.11.2)
in `reference/tests/test_functions.py`.

Root cause (confirmed with a fix, see below): `erf_impl` follows
Boost's `erf_imp` for 53-bit precision, which writes each interval's result as
exp(-z^2)/z * (Y + R(z)) with R a minimax rational fit made against Y **as a
float**: Boost declares e.g. `Y = 0.3440242112F`, which is
0.34402421116828918... statrs uses the decimal double `0.3440242112`, larger by
3.17e-11. At x = 0.5, Y + R is about 0.3078, so the predicted relative error is
3.17e-11 / 0.3078 = 1.03e-10, matching the measured 1.03e-10. All 13 interval
constants in `erf_impl` (`0.3440242112` to `0.5641584396`) differ from their
float values by 5e-12 to 4.3e-11 in the same way. Using the float-rounded values
(for example `0.3440242112_f32 as f64`) would likely fix it.

Confirmed (2026-10-14, while diagnosing [finding 009](009-statrs-normal.md)):
a copy of statrs 0.19.1 with all 13 constants replaced by `<constant>_f32 as f64`
passes all 281 points of `erfc.jsonl` (worst 3.8e-16, tolerance 1e-13).
statrs's `Normal::cdf` and `Normal::sf` inherit this error (finding 009).
