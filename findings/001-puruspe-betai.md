# 001: puruspe::betai wrong for large a and b, and inaccurate everywhere

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silently wrong answers: 0 instead of 0.5, 3% off), plus a
  low-severity accuracy loss at most other points
- **Status:** approved by owner (2026-10-14; extended scope re-confirmed the
  same day). Root cause diagnosed and a fix prepared (milestone 3); submitting
  the fix upstream is approved, not yet done.
- **Affected inputs:** the exact lists are `BETAI` and `BETAI_LARGE` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md); tolerance 1e-12.
  - `reference/data/betainc.jsonl`: 336 of 639 points.
  - `reference/data/betainc_large.jsonl` (a, b > 3000, added in milestone 3):
    **all** 144 points.
  1. **Wrong answers** (149 points).
     - **Every** input with both a and b above 3000 that is not exactly 0 or 1
       gives 0 or -0, including I_x(a, a) at x = 1/2 (true: 0.5). On the grid:
       I_x(1e4, 1e4) for x in {0.5, 0.7, 0.9, 0.99}, and all 144 points of
       `betainc_large` (a, b in {3001, 5000, 2e4, 1e5}, x from 10 standard
       deviations below the mean to 10 above).
     - I_0.9(1e4, 1e3) is 3.1% low.
  2. **Accuracy loss** (331 points of `betainc`): relative errors from 1e-12 up
     to 5.7e-10, at every parameter size on the grid, including a = b = 0.1.
     By max(a, b): 0.1 to 10: up to 3e-11; 100: up to 6e-11; 1000 and 1e4: up
     to 5.7e-10. SciPy's worst error on the same grid is 1.6e-13 (1.5e-13 on
     `betainc_large`).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // I_{1/2}(a, a) = 1/2 exactly, by the symmetry I_x(a, b) = 1 - I_{1-x}(b, a).
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.5)); // 0.0, should be 0.5
    println!("{:?}", puruspe::betai(3001.0, 3001.0, 0.5)); // 0.0, should be 0.5
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.7)); // -0.0, should be 1.0 (to double precision)
    println!("{:?}", puruspe::betai(1e4, 1e3, 0.9)); // 0.000596..., should be 0.000615...
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
4. **Accuracy loss: the prefactor.** I = exp(lnΓ(a+b) - lnΓ(a) - lnΓ(b) +
   a ln x + b ln(1-x)) times a fraction or an integral. puruspe's `ln_gamma` is
   a 7-term Lanczos fit (g = 5) whose absolute error is 0.4e-11 to 2e-11
   for arguments up to 10, 1e-10 at 100 and about 2e-10 from 1e3 to 2e4
   (measured against glibc's `lgamma`); that error passes directly into the
   result's relative error. Replacing `ln_gamma` by glibc's `lgamma` in the
   fixed code brings the failures on `betainc` from 327 to 30 (worst
   1.05e-10) and on `betainc_large` from 112 to 100 (worst 3.4e-10). The
   rest is the cancellation the tolerance reason describes: at a = 1e5 the logarithms are
   of order 1e6, whose rounding alone costs about 2e-10. Removing it needs a
   different prefactor (e.g. Boost's Lanczos-ratio `ibeta_power_terms`), which
   is beyond a minimal fix. `gammp` ([002](002-puruspe-gammp.md)) shares
   `ln_gamma` and shows the same floor.

## Fix

[fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch](../fixes/puruspe/0001-Fix-betai-for-large-a-and-b.patch)
fixes bugs 1 to 3 (three lines) and adds a regression test to puruspe's
`tests/beta_test.rs` that fails before and passes after. It applies to puruspe
0.4.4 and to upstream `main` (commit `f3eb94b`, same file). Not submitted:
that needs the owner's approval.

Verified here as `puruspe_fix::betai` (crate `crates/puruspe-fix`, a verbatim
copy of the patched items, kept identical by `scripts/check_fixes.py`); report
[reports/puruspe-fix-0.1.0.md](../reports/puruspe-fix-0.1.0.md):

| Dataset | puruspe 0.4.4: failing points, worst | with the fix: failing points, worst |
|---|---|---|
| `betainc` (639) | 336, 1.0 (wrong answers) | 327, 5.15e-10 |
| `betainc_large` (144) | 144, 1.0 (wrong answers) | 112, 4.1e-10 |

No wrong answers remain; every remaining failure is bug 4 and is listed in
`FIX_BETAI` and `FIX_BETAI_LARGE` in `crates/check-puruspe/src/expected.rs`.
**Trade-off of fix 3.** Points with small b and x close to 1 now take the
`1 - (...)` branch, which loses accuracy when the result is small. On
`betainc`, 41 points get worse by more than 1.5x (27 get better); the worst on
that grid is (100, 0.1, 0.99), from 1.3e-11 to 5.15e-10.

**Random spot-check (2026-10-14, before submission): the trade-off is larger
than the grid showed.** 5,000 random points (a, b log-uniform in 1e-2..1e5;
x uniform, near the mean, or within 1e-6 of 0 or 1), original and fixed code
compared with SciPy, every disagreement re-checked with mpmath at 50 and 100
digits:

- Wrong answers (relative error > 1e-6): 511 with puruspe 0.4.4, 2 with the
  fix. Both remaining ones are the prefactor floor (bug 4) at extreme
  parameters, (30389, 8425, 0.699) and (58131, 25112, 0.642), error 1.2e-6 and
  1.1e-6 on values near 1e-302 and 1e-256 that 0.4.4 returned as 0. Two more
  flagged points were SciPy's errors (it returns 0 or is 3e-6 off near 1e-292;
  puruspe is right).
- **Regression from fix 3: 38 points get more than 10x worse, up to 1.8e-8**
  (vs 2.4e-10 before), all with b < 0.06 and the result I between 0.003 and
  0.02. Worst, against mpmath:

  | (a, b, x) | 0.4.4 error | fixed error |
  |---|---|---|
  | (73256.2, 0.01866, 0.9999945374) | 2.4e-10 | 1.8e-8 |
  | (58386.6, 0.01965, 0.9999866779) | 1.1e-10 | 1.5e-8 |
  | (97698.9, 0.0166, 0.9999973151) | 1.8e-10 | 1.0e-8 |
  | (7457.77, 0.02365, 0.9998657733) | 3.2e-11 | 5.9e-9 |
  | (15.766, 0.01229, 0.9565398173) | 1.8e-11 | 3.7e-9 |

  Cause: for tiny b the switch point (a+1)/(a+b+2) is close to 1, so these x
  use 1 - I_{1-x}(b, a) with I_{1-x}(b, a) around 0.99; the subtraction
  magnifies the prefactor's error about 100-fold. The grid had no b below 0.1,
  so CI did not see it.

**Consequence:** the patch is **not submitted**. Its description understates
the trade-off, and fix 3 makes some previously accurate results worse. Next:
add small-b points to the reference data so CI covers this region, then
revise fix 3 (for example, choose the branch whose result is the smaller of I
and 1 - I, as Boost does) and re-verify, including a random spot-check.

## Notes

The true values rely on the identity I_{1/2}(a, a) = 1/2 where it applies, and
otherwise on the generator's positive-term series (DLMF 8.17.8), which is
tested against `mpmath.betainc` and the exact binomial sum DLMF 8.17.5
(`reference/tests/test_functions.py`). SciPy agrees with mpmath at every point.
