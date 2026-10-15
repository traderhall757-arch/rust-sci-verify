# 030: statrs ln_gamma is finite at the poles and loses all digits next to 1 and 2

- **Crate / version:** statrs 0.19.1 (`statrs::function::gamma::ln_gamma`)
- **Severity:** high (silently wrong answers: `ln_gamma(-2.0)` is 36.4 instead
  of +inf, `ln_gamma(1.0)` is -4.4e-16 instead of 0, and
  `ln_gamma(0.999999999999999)` is 0 instead of 5.8e-16)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 32 of the 877 points of `reference/data/gammaln.jsonl`
  (tolerance 2e-13) where Gamma(x) > 0 or x is a pole. The exact list is
  `LN_GAMMA` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  1. **Poles** (19 of the 36 pole points, the negative integers -2 to -30
     where it happens): finite values from -41.4 to 36.4 instead of +inf.
     At the other poles it returns NaN or an infinity, which the conventions
     accept.
  2. **Next to 1 and 2** (13 points, |x - 1| or |x - 2| from 1e-15 to 1e-5,
     and x = 1 itself): an absolute error of up to about 4e-15, so the
     relative error grows as ln Gamma(x) goes to 0: 3.2e-10 at 2.00001,
     1.5e-6 at 1.9999999999, 7.5 at 1.999999999999999, 100% (0 returned) at
     0.999999999999999.

statrs documents `ln_gamma` as "the logarithm of the gamma function", so
for x < 0 where Gamma(x) < 0 (an odd floor, e.g. -0.5) its NaN is its
documented answer; those 79 points of the dataset (which holds log|Gamma|)
are not run. Gamma is positive elsewhere, and there it passes every point
away from 1, 2 and the poles.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::function::gamma::ln_gamma;

fn main() {
    println!("{:?}", ln_gamma(-2.0)); // a pole: should be +inf
    println!("{:?}", ln_gamma(-17.0)); // a pole: should be +inf
    println!("{:?}", ln_gamma(1.0)); // should be 0.0
    println!("{:?}", ln_gamma(0.999999999999999)); // should be 5.767543092132013e-16
    println!("{:?}", ln_gamma(1.9999999999)); // should be -4.2278437004755317e-11
    println!("{:?}", ln_gamma(2.00001)); // should be 4.22787559764836e-6
    println!("{:?}", ln_gamma(-0.5)); // Gamma(-0.5) < 0: NaN, as documented (ln of Gamma)
}
```

This is the program that was run (release and debug mode give the same
output; it also printed the digamma lines of finding 031).

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| -2 | 36.397144537438486 | +inf (pole) | +inf | pole |
| -17 | 1.792616759836676 | +inf (pole) | +inf | pole |
| 1 | -4.440892098500626e-16 | 0 | 0 | 2.0e292 (relative to 2.2e-308) |
| 0.999999999999999 | 0.0 | 5.767543092132013e-16 | 5.551115123125783e-16 | 1 |
| 1.9999999999 | -4.227906913456536e-11 | -4.2278437004755317e-11 | -4.2278403000750586e-11 | 1.5e-6 |
| 2.00001 | 4.2278755962854575e-6 | 4.22787559764836e-6 | 4.227875597648359e-6 | 3.2e-10 |

Second opinion: SciPy (`scipy.special.gammaln`) agrees with mpmath at the
poles, at 1 and at 2.00001; next to 1 and 2 it is off too (3.8% at
0.999999999999999, 8.0e-7 at 1.9999999999), so there the second route is
the series ln Gamma(1 + z) = -gamma z + sum (-1)^k zeta(k) z^k / k
(DLMF 5.7.3) and ln Gamma(2 + z) = ln Gamma(1 + z) + ln(1 + z), with z
exact, tested against the dataset's route in
`reference/tests/test_functions.py` (from the `gammaln` dataset's PR). At
the poles the true value needs no computation: Gamma has a pole at every
non-positive integer (DLMF 5.2.1), so |Gamma| and its logarithm are +inf.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the `gammaln` route is tested against
   log(gamma(x)), the reflection formula, the series next to 1 and 2 and
   Stirling's series (`reference/tests/test_functions.py`).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 4,000 random x (Python `random.Random(20261019)`): 25%
   x = 1 or 2 +- 10^u with u uniform in [-16, -1]; 10% negative integers
   -1..-1e6; 35% 10^u with u uniform in [-300, 300]; 30% -10^u with u
   uniform in [-300, 15]; 6 significant digits except the first group;
   reference `mpmath.loggamma` (real part) at 50 and 100 digits. 1,165 with
   Gamma(x) < 0 are left out, as above. Results: 215 of the 435 poles give a
   finite value. Next to 1 and 2, by |x - n|: every x below 1e-5 fails
   (739 of 739), 122 of 126 from 1e-5 to 1e-3, 48 of 63 from 1e-3 to 1e-2,
   8 of 89 from 1e-2 to 0.5 (worst relative error 13, where ln Gamma is
   nearly 0). All 1,366 positive
   x away from 1 and 2 and all 9 negative x with Gamma(x) > 0 pass.
5. **Root cause** (from the source, `src/function/gamma.rs`): both
   branches evaluate the Lanczos approximation as a sum of logarithms,
   `s.ln() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / E).ln()`
   for x >= 0.5, whose terms are of size 1 to 10 and cancel to
   ln Gamma(x) = 0 at x = 1 and 2: each is rounded, so the sum carries an
   absolute error of a few 1e-16 that becomes relative error as the result
   goes to 0. Careful implementations (fdlibm, Boost, SciPy's Cephes) switch
   to a rational approximation in x - 1 or x - 2 there. For x < 0.5 the
   reflection uses `(PI * x).sin().ln()`; at a negative integer x, PI * x is
   rounded, so its sine is a tiny nonzero number (2.4e-16 at x = -2) instead
   of 0, and the result is finite instead of +inf; where the rounded sine is
   negative the logarithm gives NaN (accepted at a pole). A test for a
   non-positive integer x before the reflection would return +inf.

## Notes

The poles and the values next to 1 and 2 matter in practice: ln Gamma(1) =
0 is ln(0!), which comes out as -4.4e-16 (ln Gamma(2) = ln(1!) is exactly
0). statrs's distributions call `ln_gamma` in 17 places (`Gamma`, `Beta`,
`StudentsT`, `Chi`, `InverseGamma`, `Dirichlet` and others), mostly for
densities and entropies.
