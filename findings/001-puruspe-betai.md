# 001: puruspe::betai wrong for large a and b and near its switch point, and inaccurate everywhere

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silently wrong answers: 0 instead of 0.5, 1.0 instead of
  3.4e-14, 3% off), plus a low-severity accuracy loss at most other points
- **Status:** approved by owner (2026-10-14; extended scope re-confirmed the
  same day). Root cause diagnosed and a fix prepared (milestone 3). A random
  spot-check before submission showed that bug 3 also gives wrong answers at
  ordinary parameters (a, b <= 3000), and that the first fix for it lost
  accuracy for small b; the fix was revised (milestone 3 follow-up). Wider
  scope confirmed by owner (2026-10-14). Not submitted: nothing goes upstream
  until the project is finished and tested (owner's decision).
- **Affected inputs:** the exact lists are `BETAI`, `BETAI_LARGE` and `BETAI_SWITCH` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md); tolerance 1e-12.
  - `reference/data/betainc.jsonl`: 336 of 639 points.
  - `reference/data/betainc_large.jsonl` (a, b > 3000, added in milestone 3):
    **all** 144 points.
  - `reference/data/betainc_switch.jsonl` (small b with x near 1, and points
    between two candidate switch points; added in the milestone 3 follow-up):
    170 of 186 points.
  1. **Wrong answers** (149 points of `betainc` and `betainc_large`, 66 of
     `betainc_switch`).
     - **Every** input with both a and b above 3000 that is not exactly 0 or 1
       gives 0 or -0, including I_x(a, a) at x = 1/2 (true: 0.5). On the grid:
       I_x(1e4, 1e4) for x in {0.5, 0.7, 0.9, 0.99}, and all 144 points of
       `betainc_large` (a, b in {3001, 5000, 2e4, 1e5}, x from 10 standard
       deviations below the mean to 10 above).
     - I_0.9(1e4, 1e3) is 3.1% low.
     - Continued-fraction branch (a or b <= 3000), on the wrong side of the
       mean (bug 3): results near 1 where the true value is tiny, e.g.
       I_0.41683(1000, 1000) = 1.0 (true: 3.4e-14) and I_0.6(1000, 500) = 1.05
       (true: 4.8e-8); and, for b < 2 with x close to 1, results far too low,
       e.g. I_x(1, 0.5) = 0.197 at x = 0.9999999999 (true: 0.99999, from the
       closed form 1 - sqrt(1 - x)) and I_0.99999999(50, 0.2) 2% low. On
       `betainc_switch`: 66 points off by more than 1e-6 (60 with small b,
       6 between the switch points). The `betainc` grid has only one such
       point, (1e4, 1e3, 0.9) above.
  2. **Accuracy loss** (331 points of `betainc`): relative errors from 1e-12 up
     to 5.7e-10, at every parameter size on the grid, including a = b = 0.1.
     By max(a, b): 0.1 to 10: up to 3e-11; 100: up to 6e-11; 1000 and 1e4: up
     to 5.7e-10. SciPy's worst error on the same grid is 1.6e-13 (1.5e-13 on
     `betainc_large`, 5.9e-14 on `betainc_switch`). On `betainc_switch`, 104
     more points have errors between 1e-12 and 1e-6.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // I_{1/2}(a, a) = 1/2 exactly, by the symmetry I_x(a, b) = 1 - I_{1-x}(b, a).
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.5)); // 0.0, should be 0.5
    println!("{:?}", puruspe::betai(3001.0, 3001.0, 0.5)); // 0.0, should be 0.5
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.7)); // -0.0, should be 1.0 (to double precision)
    println!("{:?}", puruspe::betai(1e4, 1e3, 0.9)); // 0.000596..., should be 0.000615...
    println!("{:?}", puruspe::betai(1e3, 500.0, 0.6)); // 1.0536151685808233, should be 4.796...e-8
    // I_x(1, b) = 1 - (1 - x)^b, so this should be 1 - sqrt(1e-10) = 0.99999 (to 5e-12).
    println!("{:?}", puruspe::betai(1.0, 0.5, 0.9999999999)); // 0.19737493856728083
    println!("{:?}", puruspe::betai(1e3, 1e3, 0.5)); // 0.5000000000864451, should be 0.5
    println!("{:?}", puruspe::betai(0.1, 0.1, 0.5)); // 0.5000000000082913, should be 0.5
}
```

## Expected

| Input (a, b, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (1e4, 1e4, 0.5) | 0.0 | 0.5 | 0.4999999999999989 | 1.0 |
| (3001, 3001, 0.5) | 0.0 | 0.5 | 0.49999999999999967 | 1.0 |
| (1e4, 1e4, 0.7) | -0.0 | 1.0 | 1.0 | 1.0 |
| (1e4, 1e3, 0.9) | 0.0005960975545106217 | 0.000615295567297501 | 0.0006152955672975032 | 3.1e-2 |
| (1000, 1000, 0.41683333333333333) | 1.0000000000000344 | 3.444853881009224e-14 | 3.444853881009216e-14 | 2.9e13 |
| (1, 0.5, 0.9999999999) | 0.19737493856728083 | 0.9999899999995863 | 0.9999899999995863 | 0.80 |
| (1e3, 100, 0.9) | 0.14747917634280627 | 0.14747917625807305 | 0.1474791762580732 | 5.7e-10 |
| (1e3, 1e3, 0.5) | 0.5000000000864451 | 0.5 | 0.4999999999999992 | 1.7e-10 |
| (0.1, 0.1, 0.5) | 0.5000000000082913 | 0.5 | 0.4999999999999999 | 1.7e-11 |
| (0.1, 0.1, 0.01) | 0.32030825038093796 | 0.3203082503756263 | 0.3203082503756264 | 1.7e-11 |

## Root cause (milestone 3)

Three independent bugs cause the wrong answers; a fourth, separate issue
causes the accuracy loss. Each was confirmed by applying the changes one at a
time to a copy of the code and rerunning the reference data.

1. **Missing parentheses in the quadrature integrand** (`betaiapprox`, used
   when a and b both exceed `SWITCH` = 3000). puruspe computes
   `b1 * (1 - t).ln() - lnmuc` where Numerical Recipes (3rd ed., section 6.4)
   has `b1 * ((1 - t).ln() - lnmuc)`. The exponent is too small by
   (b1 - 1) ln(1 / (1 - mu)), about 6900 at a = b = 1e4 and still about 2080
   at a = b = 3001, so every term underflows to 0 and the integral is 0. This
   alone explains the 0 results.
2. **Branch chosen by the sign of an underflowed value.** `betaiapprox` returns
   `1 - ans` if `ans > 0` and `-ans` otherwise; the sign of `ans` is meant to
   say whether x is above the mean. When the tail integral underflows to +0
   (always, given bug 1; and legitimately far in the upper tail, e.g.
   I_0.7(1e4, 1e4), 57 standard deviations out), the result is -0 instead
   of 1. Choosing the branch by `x > mu` fixes it. Numerical Recipes has the
   same flaw.
3. **Wrong continued-fraction switch point** (`betai`). puruspe evaluates the
   continued fraction for I_x(a, b) when `x < (a + 1) / (a + 2b)`, and the one
   for 1 - I_(1-x)(b, a) otherwise. Numerical Recipes uses
   `(a + 1) / (a + b + 2)`, below which the fraction converges rapidly. At
   (1e4, 1e3, 0.9) puruspe's switch is 0.833, so it evaluates the second form,
   whose fraction at x' = 1 - 0.9 lies outside its good region (x' < 0.091).
   There the evaluation is numerically unstable: the fraction's value in
   double precision is 1.5e-6 off at x' = 0.1 and 1.9e-5 off at
   x' = 1 - 0.9, two ulps away (mpmath at 50 digits as the reference). The
   1.9e-5 error in the complement is 3.1% of the small result. With the
   Numerical Recipes switch the direct fraction converges in 27 terms and
   the result is within 1.8e-10.

   The random spot-check (below) showed this is one side of a wider
   problem. Both points are close to the mean a / (a + b), and each fraction
   goes wrong on the far side of it. puruspe's switch lies below the mean
   when 1/a + 1/b < 1 (for example a, b > 2), so the complement's fraction
   is used too early: up to a factor of 1e37 off at (2900, 2900, 0.4167),
   1.0 instead of 1.2e-37. Otherwise (for example b < 1) it lies above the
   mean, and above 1 when b < 1/2, so the direct fraction is used all the
   way to x = 1: I_x(1, 0.5) at x = 0.9999999999 is 0.197 instead of
   0.99999. **Fix 3 (revised):** switch at the mean, `x < a / (a + b)`, as
   Boost.Math's `ibeta_imp` does for a, b >= 1.
4. **Accuracy loss: the prefactor.** I = exp(lnΓ(a+b) - lnΓ(a) - lnΓ(b) +
   a ln x + b ln(1-x)) times a fraction or an integral. puruspe's `ln_gamma` is
   a 7-term Lanczos fit (g = 5) whose absolute error is 0.4e-11 to 2e-11
   for arguments up to 10, 1e-10 at 100 and about 2e-10 from 1e3 to 2e4
   (measured against glibc's `lgamma`); that error passes directly into the
   result's relative error, and is magnified by (1 - I)/I when a small I is
   computed as 1 - (complement). Replacing `ln_gamma` by glibc's `lgamma` in
   the first version of the fix brings the failures on `betainc` from 327 to 30 (worst
   1.05e-10) and on `betainc_large` from 112 to 100 (worst 3.4e-10). The
   rest is the cancellation the tolerance reason describes: at a = 1e5 the logarithms are
   of order 1e6, whose rounding alone costs about 2e-10. Removing it needs a
   different prefactor (e.g. Boost's Lanczos-ratio `ibeta_power_terms`), which
   is beyond a minimal fix. `gammp` ([002](002-puruspe-gammp.md)) shares
   `ln_gamma` and shows the same floor.

## Fix

[fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch](../fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch)
fixes bugs 1 to 3 (three lines) and adds a regression test to puruspe's
`tests/beta_test.rs` that fails before and passes after (and also fails
with the first version of fix 3). It applies to puruspe 0.4.4 and to
upstream `main` (commit `f3eb94b`, same file). Not submitted: held until the
project is finished and tested, then needs the owner's approval.

Verified here as `puruspe_fix::betai` (crate `crates/puruspe-fix`, a verbatim
copy of the patched items, kept identical by `scripts/check_fixes.py`); report
[reports/puruspe-fix-0.1.0.md](../reports/puruspe-fix-0.1.0.md):

| Dataset | puruspe 0.4.4: failing points, worst | with the fix: failing points, worst |
|---|---|---|
| `betainc` (639) | 336, 1.0 (wrong answers) | 321, 1.8e-10 |
| `betainc_large` (144) | 144, 1.0 (wrong answers) | 112, 4.1e-10 |
| `betainc_switch` (186) | 170, 8.6e36 (wrong answers) | 144, 3.0e-9 |

No wrong answers remain; every remaining failure is bug 4 and is listed in
`FIX_BETAI`, `FIX_BETAI_LARGE` and `FIX_BETAI_SWITCH` in
`crates/check-puruspe/src/expected.rs`. CI also bounds the worst error per
dataset (2.5e-10, 5e-10, 3.5e-9). The worst point, (1e5, 0.01,
0.99999999), has I = 0.061 computed as 1 - 0.939: the prefactor's 2e-10
error at a = 1e5, magnified 15-fold. puruspe 0.4.4 is 1.7e-6 off there; the
direct fraction is not usable that far above the mean.

### History: the first fix 3 and the random spot-checks

The first version of fix 3 (milestone 3) used Numerical Recipes' switch
`(a + 1) / (a + b + 2)`. It removed the wrong answers, but it lies below
the mean when b < a, so for small b and x close to 1 a small I came from
1 - (a value close to 1). A random spot-check before submission
(2026-10-14; 5,000 points) found it up to 1.8e-8 off there, against
2.4e-10 for 0.4.4; worst points, with mpmath:

| (a, b, x) | 0.4.4 error | first fix | revised fix |
|---|---|---|---|
| (73256.2, 0.01866, 0.9999945374) | 2.4e-10 | 1.8e-8 | 1.0e-11 |
| (58386.6, 0.01965, 0.9999866779) | 1.1e-10 | 1.5e-8 | 2.3e-10 |
| (97698.9, 0.0166, 0.9999973151) | 1.8e-10 | 1.0e-8 | 2.1e-10 |
| (7457.77, 0.02365, 0.9998657733) | 3.2e-11 | 5.9e-9 | 2.8e-11 |
| (15.766, 0.01229, 0.9565398173) | 1.8e-11 | 3.7e-9 | 1.8e-11 |

The first two columns are from the spot-check's full-precision inputs; the
revised fix's errors are at the inputs as quoted (now points of
`betainc_switch`), where 0.4.4 has the same errors. The grid had no b
below 0.1, so CI did not see it.

The follow-up added `betainc_switch` and compared four switch rules on
20,000 new random points (a, b log-uniform in 1e-2..1e5; x uniform, near
the mean, or within 1e-6 of 0 or 1; seed 2). Each was compared with SciPy,
and every difference above 1e-9 where a or b <= 3000 re-checked with mpmath
at 50 digits (550 points), as were the 7 above 1e-9 where both exceed 3000.
Of the 18,991 points with a or b <= 3000, 3,342 have results below 1e-300,
where SciPy itself is unreliable; the table counts the other 15,649:

| Switch rule | points off by > 1e-6 | > 1e-9 | worst |
|---|---|---|---|
| 0.4.4: (a + 1) / (a + 2b) | 414 | 783 | 5e113 |
| first fix: (a + 1) / (a + b + 2) | 0 | 127 | 1.1e-8 |
| larger of the two | 146 | 304 | 0.47 |
| **revised fix: a / (a + b)** | **0** | **2** | **1.1e-9** |

With the revised fix, 44 points are more than 10 times less accurate than
in 0.4.4, all with b < 0.12 and x > 0.993, worst 7.2e-10 (was 4.3e-11):
past the mean, a small I is computed as 1 - (complement), as above. Where
a and b both exceed 3000 (1,009 points; quadrature, untouched by fix 3), 7
points are off by more than 1e-9, up to 4.8e-6, all with results below
1e-115 where 0.4.4 returned 0: the prefactor's error grows with the size of
its terms there (bug 4).

## Notes

The true values rely on the identity I_{1/2}(a, a) = 1/2 where it applies, and
otherwise on the generator's positive-term series (DLMF 8.17.8), which is
tested against `mpmath.betainc` and the exact binomial sum DLMF 8.17.5
(`reference/tests/test_functions.py`). SciPy agrees with mpmath at every point.
