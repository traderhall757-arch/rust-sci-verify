# 033: puruspe ln_gamma is NaN or wrong for negative x, inf for tiny x, and loses digits near 1 and 2

- **Crate / version:** puruspe 0.4.4 (`puruspe::ln_gamma`)
- **Severity:** high (silently wrong answers: `ln_gamma(-3.75)` is 1.8e-4 off,
  `ln_gamma(1.0)` is -5.2e-12 instead of 0, `ln_gamma(0.999999999999999)` has
  the wrong sign; plus NaN and +inf where the value is finite)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 666 of the 877 points of `reference/data/gammaln.jsonl`
  (tolerance 2e-13) where Gamma(x) > 0 or x is a pole. The exact list is
  `LN_GAMMA` in `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Four groups, by
  cause (section "Root cause"):
  1. **Negative x** (all 46 points with x < 0 and Gamma(x) > 0): NaN for
     x < -4.5 (40 points, -5.25 down to -1e15 - 0.5), and finite values off
     by 1.4e-9 to 1.8e-4 on (-4, -3) and (-2, -1) (6 points).
  2. **Tiny x** (30 points): +inf for x = 1e-300 to 1e-20 (29 points; the
     true values are 46 to 691), and 3.6e-9 off at x = 1e-10.
  3. **Next to 1 and 2** (14 points, |x - 1| or |x - 2| from 1e-15 to 1e-5,
     and x = 1 and 2 themselves): `ln_gamma(1)` = -5.2e-12 and
     `ln_gamma(2)` = -1.0e-11 instead of 0, so the relative error grows as
     ln Gamma(x) goes to 0: 9.1e-7 at 0.99999, 9% at 0.9999999999, and the
     wrong sign at 0.999999999999999.
  4. **Accuracy floor** (576 points, x = 0.25 to 155.75 away from 1 and 2):
     an absolute error of up to 1.3e-10, a relative error of 2.0e-13 to
     1.1e-10. Above x = 156 ln Gamma exceeds 600 and the same absolute error
     is below the tolerance: none of the 135 grid points from 156 to the
     largest double fails (130 pass; at the 5 where ln Gamma overflows it
     returns +inf, which the conventions accept). 45 smaller x pass, where
     the error happens to be small.

At the 31 poles it returns NaN or an infinity, which the conventions accept.
puruspe documents `ln_gamma` as "the natural logarithm of the gamma
function", so where Gamma(x) < 0 (x < 0 with an odd floor, e.g. -0.5) the
dataset's log|Gamma| is not its contract; those 79 points are not run.

`ln_gamma` is also the source of the accuracy floor of `betai`
([001](001-puruspe-betai.md)) and `gammp` ([002](002-puruspe-gammp.md)).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // ln_gamma (finding 033)
    println!("{:?}", puruspe::ln_gamma(1.0)); // should be 0.0
    println!("{:?}", puruspe::ln_gamma(2.0)); // should be 0.0
    println!("{:?}", puruspe::ln_gamma(0.999999999999999)); // should be 5.767543092132013e-16
    println!("{:?}", puruspe::ln_gamma(0.5)); // should be 0.5723649429247001
    println!("{:?}", puruspe::ln_gamma(155.0)); // should be 625.128656730891
    println!("{:?}", puruspe::ln_gamma(1e-10)); // should be 23.025850929882736
    println!("{:?}", puruspe::ln_gamma(1e-20)); // should be 46.051701859880914
    println!("{:?}", puruspe::ln_gamma(-1.5)); // should be 0.860047015376481
    println!("{:?}", puruspe::ln_gamma(-5.5)); // should be -4.517832174007741
    // gammp for large a (finding 002)
    println!("{:?}", puruspe::gammp(1e8, 99990000.0)); // should be 0.15865525352814383
    println!("{:?}", puruspe::gammp(1e8, 9e7)); // should be 0.0
    println!("{:?}", puruspe::gammp(1e9, 1e9)); // should be 0.5000042052208701
    println!("{:?}", puruspe::gammq(100.0, 50.0)); // should be 0.9999999996799934
}
```

This is the program that was run (release and debug mode give the same
output; the last four lines belong to finding 002).

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| -5.5 | NaN | -4.517832174007741 | -4.517832174007742 | inf |
| -1.5 | 0.8600470105518334 | 0.860047015376481 | 0.860047015376481 | 5.6e-9 |
| 1e-20 | inf | 46.051701859880914 | 46.051701859880914 | inf |
| 1e-10 | 23.0258508471237 | 23.025850929882736 | 23.025850929882733 | 3.6e-9 |
| 1 | -5.232259070453438e-12 | 0 | 0 | 2.4e296 (relative to 2.2e-308) |
| 0.999999999999999 | -5.232259070453438e-12 | 5.767543092132013e-16 | 5.551115123125783e-16 | 9.1e3 |
| 0.5 | 0.5723649429204194 | 0.5723649429247001 | 0.5723649429247 | 7.5e-12 |
| 155 | 625.1286567310166 | 625.128656730891 | 625.1286567308911 | 2.0e-13 |

Second opinion: SciPy (`scipy.special.gammaln`) agrees with mpmath within
1.8e-15 at every row but 0.999999999999999, where it is off by 3.8% too.
There, and at the 13 other points next to 1 and 2, the second route is the
series ln Gamma(1 + z) = -gamma z + sum (-1)^k zeta(k) z^k / k (DLMF 5.7.3)
and ln Gamma(2 + z) = ln Gamma(1 + z) + ln(1 + z) with z exact, tested
against the dataset's route in `reference/tests/test_functions.py` (as for
finding 030).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `puruspe = "=0.4.4"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the `gammaln` route is tested against
   log(gamma(x)), the reflection formula, the series next to 1 and 2 and
   Stirling's series (`reference/tests/test_functions.py`).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 4,000 random x (Python `random.Random(20261020)`): 25%
   x = 1 or 2 +- 10^u with u uniform in [-16, -1]; 10% negative integers
   -1..-1e6; 30% 10^u with u uniform in [-300, 300]; 15% uniform in
   (0, 200); 20% -10^u with u uniform in [-300, 15]; 6 significant digits
   except the first group; reference `mpmath.loggamma` (real part) at 50 and
   100 digits, which agree at every point. 713 with Gamma(x) < 0 are left
   out, as above. Results, against the tolerance 2e-13:
   - poles: 0 of 454 give a finite value;
   - negative x with Gamma(x) > 0: 9 of 9 fail (7 NaN below -4.5; 2 finite,
     worst 2.9e-4 at -3.79884);
   - x below 2^-53: all 572 give +inf (the largest of them was 3.6e-17;
     the smallest x with a finite result was 4.7e-16); from 2^-53 to 1e-2:
     24 of 24 fail, worst 1.5e-3;
   - |x - 1| or |x - 2| < 1e-2: 984 of 984 fail;
   - other x from 1e-2 to 200: 530 of 685 fail, worst 2.5e-9 (at 1.98997,
     where ln Gamma is 4e-3); with |x - 1| and |x - 2| at least 0.5, 447 of
     602 fail, worst 3.4e-11; the largest absolute error is 1.4e-10 and the
     largest failing x is 155.889;
   - x above 200: 0 of 559 fail.
5. **Root cause** (from the source, `src/gamma.rs`): `ln_gamma` evaluates a
   Lanczos approximation with g = 5 and 7 coefficients (`LG5N7`),
   `ln sqrt(2 pi) + ln s - base + (z + 0.5) ln base` with z = x - 1,
   base = z + 5.5 and s = c0 + sum c_i / (z + i), for every x, with no
   special cases. Numerical Recipes (2nd ed., section 6.1) quotes an error
   below 2e-10 for its approximation with the same g and length.
   - **Group 4**: the approximation's own error. A relative error e in Gamma
     is an absolute error e in ln Gamma, so it is above the tolerance until
     ln Gamma(x) exceeds e / 2e-13, about 600 (x = 156), which matches the
     grid and the scan.
   - **Group 3**: the same error does not vanish at the zeros of ln Gamma:
     `ln_gamma(1)` = -5.2e-12, so the relative error is unbounded next to
     1 and 2. Careful implementations (fdlibm, Boost) use a rational
     approximation in x - 1 or x - 2 there.
   - **Group 1**: there is no reflection formula. The approximation is
     evaluated outside its range: for x < -4.5, base = x + 4.5 is negative
     and `base.ln()` is NaN; on (-4, -3) and (-2, -1) it gives a finite,
     wrong value.
   - **Group 2**: z = x - 1 rounds to -1 for x below 2^-54 (5.6e-17), so
     the first term c1 / (z + 1) divides by 0 and the result is +inf. Above
     that, z + 1 differs from x by up to 2^-54, a relative error of about
     1.1e-16 / x in the dominant term (3.6e-9 at x = 1e-10). Forming the first term from x itself avoids it.

## Notes

The leads for a fix are standard: a longer Lanczos set (g = 6.02, 13 terms,
as in Boost) or Stirling's series for large x, rational approximations next
to 1 and 2, and the reflection formula ln|Gamma(x)| = ln(pi / |sin(pi x)|) -
ln Gamma(1 - x) with an exactly reduced argument for x < 0. A better
`ln_gamma` is the core of the planned accuracy milestone for `betai` and
`gammp` (roadmap item 5).
