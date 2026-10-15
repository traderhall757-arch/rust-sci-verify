# 031: statrs digamma loses all digits next to its zero and digits at large negative x

- **Crate / version:** statrs 0.19.1 (`statrs::function::gamma::digamma`)
- **Severity:** high (silently wrong answers: next to the zero
  x0 = 1.46163 the sign and size are wrong, e.g. -1.4e-15 instead of
  -9.2e-17 at the double nearest x0, and at x = -999999999999999.5 it
  returns 35.07 instead of 34.54, 1.5% off)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 9 of 949 points of `reference/data/digamma.jsonl`
  (tolerance 4e-13). The exact list is `DIGAMMA` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  1. **Next to the zero x0 = 1.4616321449683622** (5 points, x0 and
     x0 +- 1e-10, 1e-5): an absolute error of about 1.3e-15, so a relative
     error of about 1.3e-15 / |psi(x)|: 14 at x0, 1.4e-5 at x0 +- 1e-10,
     1.3e-10 at x0 +- 1e-5.
  2. **Negative x** (4 points): 1.5% off at -999999999999999.5, 6.5e-11 at
     -1000000.5, 4.9e-12 at -22.75 and 1.1e-12 at -21.75.

Unlike special's `digamma` (finding 017) it does not recurse, so it cannot
overflow the stack; it returns -inf at the non-positive integers and passes
every other point, including all positive x away from x0.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::function::gamma::digamma;

fn main() {
    println!("{:?}", digamma(1.4616321449683622)); // should be -9.241265521729427e-17
    println!("{:?}", digamma(1.4616321450683623)); // should be 9.676714013423533e-11
    println!("{:?}", digamma(-999999999999999.5)); // should be 34.538776394910684
    println!("{:?}", digamma(-1000000.5)); // should be 13.815511557963816
}
```

This is the program that was run (release and debug mode give the same
output; it ran together with finding 030's lines).

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 1.4616321449683622 (nearest x0) | -1.406535478365445e-15 | -9.241265521729427e-17 | -9.241265521729427e-17 | 14 |
| 1.4616321450683623 (x0 + 1e-10) | 9.676579600537832e-11 | 9.676714013423533e-11 | 9.676714013423533e-11 | 1.4e-5 |
| -999999999999999.5 | 35.07050770235261 | 34.538776394910684 | 34.538776394910684 | 1.5e-2 |
| -1000000.5 | 13.81551155885591 | 13.815511557963816 | 13.815511557963816 | 6.5e-11 |

Second opinion: SciPy (`scipy.special.digamma`) agrees with mpmath to the
last digit at every row.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; `mpmath.digamma` is tested against the
   recurrence, the reflection formula and values at rational points in
   `reference/tests/test_functions.py` (from the `digamma` dataset's PR,
   finding 017).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 4,000 random x (Python `random.Random(20261019)`, after
   finding 030's 4,000): 25% x0 +- 10^u with u uniform in [-16, -1]; 35%
   10^u with u uniform in [-300, 300]; 40% -10^u with u uniform in
   [-300, 15]; 6 significant digits except the first group; non-positive
   integers skipped; reference `mpmath.digamma` at 50 and 100 digits. 944
   of the 1,033 points near x0 fail (worst 14), 12 of 1,524 negative x (worst
   1e-11), none of the 1,385 positive x away from x0.
5. **Root cause** (from the source, `src/function/gamma.rs`, AS 103):
   - **Cause 1:** for x < 12 the result is sum_k -1/(x + k) plus the
     asymptotic series at x + k >= 12, terms of size up to about 2.5 that
     cancel to psi(x0) = 0; their rounding leaves an absolute error of about
     1e-15. There is no expansion around x0 as in Boost and SciPy.
   - **Cause 2:** for x < 0 it uses `digamma(1.0 - x) + PI / (-PI * x).tan()`;
     PI * x is rounded, an absolute error of up to |PI x| 2^-53 in the
     argument of tan (0.35 at x = -1e15), which shifts the cotangent.
     Confirmed: forming the reflection term from the exactly reduced
     argument r = x - round(x) (exact for these x) makes all 12 failing
     negative scan points and all 4 failing negative grid points pass (0 of
     the 125 negative grid points fail); it does not change positive x, so
     the 944 points near x0 still fail.

## Notes

The same two causes are in special's `digamma` (finding 017, cause 1) and
`trigamma` (finding 018, cause 2). Within statrs, `digamma` is called by the
`Beta`, `Gamma`, `InverseGamma`, `StudentsT`, `Chi` and `Dirichlet`
distributions (entropy and related formulas), not by any cdf or quantile.
