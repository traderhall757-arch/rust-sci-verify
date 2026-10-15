# 002: puruspe::gammp returns NaN for a >= 100, and is inaccurate elsewhere

- **Crate / version:** puruspe 0.4.4
- **Severity:** medium (NaN instead of a value; no silently wrong large errors),
  plus a low-severity accuracy loss
- **Status:** approved by owner (2026-10-14); extended to `gammainc_large`
  and root cause found 2026-10-15 (Claude, checklist complete; section
  "Extension")
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

  Also 25 of 51 points of `reference/data/gammainc_large.jsonl` (a >= 1e8,
  tolerance 1e-12; list `GAMMP_LARGE`, added 2026-10-15):
  4. **NaN** (19 points): every grid point with 0.2a <= x <= a - 1, at
     every shape from 1e8 to 1e20, including x = a - sqrt(a) (true
     P = 0.159) and the points where P rounds to 0. In general: NaN for
     every a >= 100 and 0.2a <= x <= a - 1 (section "Root cause").
  5. **Accuracy loss near the mean** (6 points, x = a, a + sqrt(a) and
     a + 3 sqrt(a) at a = 1e8 and 1e9): up to 1.0e-7 (a = 1e8) and
     1.5e-6 (a = 1e9). Off the grid it reaches 1.3e-5 (a from 1e9 to
     1e10).
  `gammq` uses the same branch: `gammq(100, 50)` is NaN as well (not
  checked against a dataset; it is in the reproduction).

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

Second opinion for the `gammainc_large` points: SciPy (`scipy.special.gammainc`)
agrees with mpmath within 2.7e-15 except at x = a - 10 sqrt(a) and
a - 30 sqrt(a) for a = 1e8 and 1e9, where it is off too (the 4 "both wrong"
points); there numerical integration of the defining integral confirms the
dataset's route (`test_gammainc_large_matches_numerical_integration` in
`reference/tests/test_functions.py`).

| Input (a, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (1e8, 99990000) | NaN | 0.15865525352814383 | 0.15865525352814386 | inf |
| (1e8, 9e7) | NaN | 0.0 | 0.0 | inf |
| (1e8, 99700000) | NaN | 1.990882289646706e-198 | 1.9861345555962874e-198 | inf |
| (1e9, 1e9) | 0.5000034596065188 | 0.5000042052208701 | 0.5000042052208701 | 1.5e-6 |

The reproduction for these is in [finding 033](033-puruspe-ln-gamma.md)
(one program, run in a fresh `cargo new` project with `puruspe = "=0.4.4"`,
release and debug mode).

## Root cause (2026-10-15)

From the source (`src/gamma.rs`, `gammpapprox`, used for a >= 100 and
x >= 0.2a), confirmed by experiment:

- **NaN** (groups 1 and 4): Numerical Recipes' `gammpapprox` integrates
  from x to an end point `xu`; below the mean (x <= a - 1) `xu` is
  max(0, min(a - 1 - 7.5 sqrt(a - 1), x - 5 sqrt(a - 1))), which is less
  than x. There the product `sum * (xu - x) * exp(..)` is negative and
  Numerical Recipes returns P = -ans. puruspe's log-sum-exp rewrite takes
  `(xu - x).ln()` of the negative number, so the result is NaN for every
  such point, whatever P is. A Python copy of the function that uses
  ln|xu - x| and P = ans has no NaN on 205 random points of that region
  (`random.Random(20261025)`: a = 10^u with u uniform in [2, 9],
  x = a + z sqrt(a) with z uniform in [-8, 0], plus the grid's NaN points;
  true values below 1e-290 left out); its worst error is 2.9e-6, the floor
  below.
- **Accuracy loss near the mean** (group 5): two causes of the same size.
  The scale `a1 * (lna1 - 1) - gln` subtracts terms of size a ln a (2e10
  at a = 1e9) whose spacing of doubles is about 4e-6, plus `ln_gamma`'s own
  error ([033](033-puruspe-ln-gamma.md)); and the 18-point Gauss-Legendre
  rule is not exact. A Python copy of `gammpapprox` (it matches the crate
  bit for bit at 56 of 60 random points with a from 1e5 to 1e10 and x from
  the mean to 5 sd above, `random.Random(20261024)`, and within 6.3e-16 at
  the other 4, where Python's and Rust's `exp`/`ln` differ) has a worst
  error of 1.3e-6; with the scale computed exactly 3.6e-6, with the
  integral computed exactly 2.6e-6, with both exact 0.
- **Accuracy loss elsewhere** (group 2): the leads below.

## Extension (2026-10-15): confirmation checklist

1. **Reproduced**: see above.
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits
   to round to the same double; the `gammainc_large` route (the positive
   series DLMF 8.7.1 and Legendre's continued fraction DLMF 8.9.2) is tested
   against numerical integration at the bulk points, including the four
   where SciPy is off.
3. **Second opinion**: see above.
4. **Off-grid**: random points with a = 10^u (6 significant digits) and
   x = a + z sqrt(a), z uniform in [-40, 40] (60%), or x = a 10^v, v
   uniform in [-3, 0.3] (40%); reference by the dataset's route at 30
   digits:
   - `random.Random(20261021)`, 1,000 points, u uniform in [2, 6]: NaN at
     all 354 points with 0.2a <= x <= a - 1 and nowhere else; with
     x > a - 1, 21 of 315 fail, worst 1.5e-10; with x < 0.2a, 45 of 331
     fail (groups 2 and 3, worst 1.7e-10 among values above 1e-306;
     results up to 3.0e-307 are returned as 0);
   - `random.Random(20261022)`, 300 points, u uniform in [2, 10]: NaN at all
     108 points of that region and nowhere else; 6 of 92 fail above it
     (worst 4.8e-10), 14 of 100 below 0.2a (worst 1.7e-10);
   - `random.Random(20261023)`, 200 points near the mean, u uniform in
     [5, 10], x = a + z sqrt(a) with z uniform in [0, 5]: 154 fail; worst
     per decade of a 3.1e-10 (1e5), 9.2e-9, 3.1e-8, 7.6e-7 and 1.3e-5
     (1e9 to 1e10).
5. **Root cause**: see above. Severity stays medium: NaN instead of a value,
   and an accuracy loss.

## Notes

The true values come from `mpmath.gammainc(a, 0, x, regularized=True)`, tested
against the finite sum Q(a, x) = e^-x sum_{k<a} x^k/k! for integer a (DLMF
8.4.10), including P(500, 100), in `reference/tests/test_functions.py`.

Leads for the accuracy loss of group 2 (the NaN is diagnosed above):

- The accuracy floor of about 1e-11 matches the one in `betai`
  ([001](001-puruspe-betai.md)); both use the same short Lanczos `ln_gamma`.
