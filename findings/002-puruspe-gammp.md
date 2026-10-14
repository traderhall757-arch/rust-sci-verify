# 002: puruspe::gammp returns NaN for a >= 100, and is inaccurate elsewhere

- **Crate / version:** puruspe 0.4.4
- **Severity:** medium (NaN instead of a value; no silently wrong large errors),
  plus a low-severity accuracy loss
- **Status:** approved by owner (2026-10-14)
- **Affected inputs:** 78 of 192 points of `reference/data/gammainc.jsonl`
  (tolerance 1e-12); the exact list is `GAMMP` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md).
  1. **NaN** (8 points), all with a >= 100 and 0.2a <= x < a:
     (100, 20), (100, 50), (200, 50), (200, 100), (500, 100), (500, 200),
     (1000, 200), (1000, 500). The true values range from 3.2e-10 down to 0
     (P(1000, 200) underflows).
  2. **Accuracy loss** (69 points): relative errors from 1e-12 up to 4.9e-11, for
     every a on the grid (0.1 to 1000). SciPy's worst error on the same grid is
     3.5e-13.
  3. **Subnormal result flushed to zero** (1 point): P(200, 2) = 2.785e-316 is
     returned as 0. SciPy also returns 0 here (verdict "both wrong"); under
     [the error metric](../docs/conventions.md#the-error-metric) that is a
     relative error of 1.3e-8. Low severity.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::gammp(500.0, 100.0)); // NaN, should be 3.8087542289264636e-178
    println!("{:?}", puruspe::gammp(100.0, 50.0)); // NaN, should be 3.200065324585125e-10
    println!("{:?}", puruspe::gammp(50.0, 50.0)); // 0.518808315446501, should be 0.5188083154720433
}
```

## Expected

| Input (a, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (100, 20) | NaN | 3.488878669689653e-37 | 3.488878669689553e-37 | inf |
| (100, 50) | NaN | 3.200065324585125e-10 | 3.2000653245851495e-10 | inf |
| (500, 100) | NaN | 3.8087542289264636e-178 | 3.8087542289263496e-178 | inf |
| (1000, 500) | NaN | 3.2982727970671e-86 | 3.298272797068251e-86 | inf |
| (50, 50) | 0.518808315446501 | 0.5188083154720433 | 0.5188083154720433 | 4.9e-11 |
| (0.1, 1) | 0.9758726562860744 | 0.9758726562736723 | 0.9758726562736721 | 1.3e-11 |
| (200, 2) | 0.0 | 2.7852507e-316 | 0.0 | 1.3e-8 |

## Notes

The true values come from `mpmath.gammainc(a, 0, x, regularized=True)`, tested
against the finite sum Q(a, x) = e^-x sum_{k<a} x^k/k! for integer a (DLMF
8.4.10), including P(500, 100), in `reference/tests/test_functions.py`.

Root cause: not yet diagnosed. Leads from reading the source:

- For a >= 100 (`ASWITCH`) and x >= 0.2a the function uses `gammpapprox`
  (Gauss-Legendre quadrature); every NaN lies in that branch.
- The accuracy floor of about 1e-11 matches the one in `betai`
  ([001](001-puruspe-betai.md)); both use the same short Lanczos `ln_gamma`.
