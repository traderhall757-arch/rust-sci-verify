# 003: statrs gamma is finite at its poles, overflows early, and has 2.5e-13 error

- **Crate / version:** statrs 0.19.1 (`statrs::function::gamma::gamma`)
- **Severity:** high (silently wrong: a large finite number where the function
  has a pole), plus medium (infinity for finite values) and low (accuracy)
- **Status:** approved by owner (2026-10-14)
- **Affected inputs:** 187 of 901 points of `reference/data/gamma.jsonl`
  (tolerance 1e-13); the exact list is `GAMMA` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  1. **Finite values at poles** (30 points): every negative integer from -1 to
     -30 returns a finite number, e.g. gamma(-1) = -2.57e16, gamma(-6) = 5.94e12,
     gamma(-18) = 0.22. Only gamma(0) returns infinity. Accepted results at a pole
     are NaN, an infinity or an error ([conventions](../docs/conventions.md)).
  2. **Premature overflow** (8 points): x = 169.75 to 171.5 returns +infinity,
     although gamma(x) is finite up to x = 171.62 (gamma(171.5) = 9.48e307).
  3. **Accuracy** (149 points): relative errors from 1e-13 up to 2.5e-13 for
     x in 10..160 and for some negative x. The documentation promises "an
     accuracy of 16 floating point digits"; SciPy's worst error on the same grid
     is 5.8e-16.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::function::gamma::gamma;

fn main() {
    println!("{:?}", gamma(-6.0)); // 5938206200927.656, should be a pole (NaN or infinity)
    println!("{:?}", gamma(-1.0)); // -2.5653050788007548e16, should be a pole
    println!("{:?}", gamma(171.5)); // inf, should be 9.4833675668248e307
    println!("{:?}", gamma(65.0)); // 1.268869321859162e89, should be 64! = 1.2688693218588417e89
}
```

## Expected

| Input x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| -1 | -2.5653050788007548e16 | pole | NaN | pole |
| -6 | 5938206200927.656 | pole | NaN | pole |
| -18 | 0.2226001720186322 | pole | NaN | pole |
| -30 | 1.0987439809047159e-18 | pole | NaN | pole |
| 169.75 | inf | 1.18337011830951e304 | 1.18337011830951e304 | inf |
| 171.5 | inf | 9.4833675668248e307 | 9.483367566824803e307 | inf |
| 65 | 1.268869321859162e89 | 1.2688693218588417e89 | 1.2688693218588415e89 | 2.5e-13 |

## Notes

Gamma has simple poles at 0, -1, -2, ... (DLMF 5.2.1). The true values come
from `mpmath.gamma`, tested against gamma(n) = (n-1)! and the half-integer
formula in `reference/tests/test_functions.py`. SciPy agrees with mpmath at
every non-pole point.

Root cause (from the source, not yet tested with a fix):

- Poles: for x < 0.5 the code uses the reflection formula
  pi / (sin(pi x) * ...). At a negative integer, `(PI * x).sin()` is not exactly
  zero because pi*x is rounded, so the result is a huge or tiny finite number
  instead of infinity. No check for non-positive integers precedes it.
- Overflow: it computes `((x - 0.5 + r)/e).powf(x - 0.5)` and multiplies by the
  Lanczos sum afterwards; the power overflows before the (small) factor brings it
  back into range.
- Accuracy: the same power, with the base rounded to a double and raised to
  x - 0.5, plausibly accounts for part of the error; not analysed further.
