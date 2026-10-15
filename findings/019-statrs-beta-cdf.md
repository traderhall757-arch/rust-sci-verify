# 019: statrs Beta cdf returns exactly 1 within 1e-15 of x = 1, and cdf and sf lose digits for large shapes

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Beta`, methods
  `cdf` and `sf` of `ContinuousCDF`; both call
  `statrs::function::beta::beta_reg`)
- **Severity:** high (a silently wrong answer: `cdf` returns exactly 1 for
  0.917 at a = 1e4, b = 0.1, x = 1 - 1e-15, 9.1% off), plus an accuracy loss
  that grows with the larger shape: up to 5.0e-10 on the grid and 1.5e-9 off
  it
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 156 of 1895 points of `reference/data/beta_cdf.jsonl`
  and 129 of 1895 points of `reference/data/beta_sf.jsonl` (tolerance
  1e-12); the exact lists are `BETA_CDF` and `BETA_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two groups, by
  cause (section "Root cause"):
  1. **Exactly 1 next to x = 1** (19 `cdf` points, none of `sf`): at
     x = 1 - 1e-15 (the double 0.999999999999999) `cdf` returns 1 for every
     a when b = 0.1 or 0.5, and for a = 1e4 when b = 1. The error is the
     true upper tail 1 - I: from 1.0e-11 (a = 1e4, b = 1) and 5.6e-9
     (a = 0.1, b = 0.5) up to 9.1e-2 (a = 1e4, b = 0.1). Every x within 1e-15
     of 1 gives exactly 1; that is wrong whenever the upper tail there is
     above about 1e-16, i.e. for b below about 1.
  2. **Large shapes** (137 `cdf` and all 129 `sf` points; every one has
     max(a, b) >= 100): errors from 1e-12 to 6e-11 across the bulk and both
     tails, growing with the larger shape (worst per value of max(a, b):
     9.9e-12 at 100, 2.8e-11 at 1000, 5.0e-10 at 1e4), largest where the
     result is formed as 1 minus the other tail: 3.2e-10 for
     `cdf(0.9998951379361433)` at a = 1e4, b = 0.1 and 5.0e-10 for its
     mirror image, `sf(0.00010486206385669913)` at a = 0.1, b = 1e4.

statrs's `StudentsT` uses the same `beta_reg` and shares both causes
([finding 010](010-statrs-students-t.md), causes 2 and 3).
`Beta::inverse_cdf` has its own, worse problems:
[finding 020](020-statrs-beta-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Beta, ContinuousCDF};

fn main() {
    let beta = |a, b| Beta::new(a, b).unwrap();
    println!("{:?}", beta(1e4, 0.1).cdf(0.999999999999999)); // should be 0.9165122714031607
    println!("{:?}", beta(0.5, 0.5).cdf(0.999999999999999)); // should be 0.9999999798763621
    println!("{:?}", beta(79372.2, 0.106553).cdf(0.9999999999999997)); // should be 0.921339739174892
    println!("{:?}", beta(1e4, 1e4).cdf(0.5)); // should be 0.5
    println!("{:?}", beta(1e4, 1e4).cdf(0.46464554479100484)); // should be 6.739002403910842e-24
    println!("{:?}", beta(1.0, 1e4).cdf(1e-100)); // should be 1.0000000000000001e-96
    println!("{:?}", beta(1e4, 0.1).cdf(0.9998951379361433)); // should be 0.022330569703958216
    println!("{:?}", beta(0.1, 1e4).sf(0.00010486206385669913)); // should be 0.022330569703971435
    println!("{:?}", beta(71839.6, 0.131114).cdf(0.9999877675)); // should be 0.039496308017432535
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| a | b | x | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|---|
| 1e4 | 0.1 | 1 - 1e-15 | `cdf` | 1.0 | 0.9165122714031607 | 0.9165122714031606 | 9.1e-2 |
| 0.5 | 0.5 | 1 - 1e-15 | `cdf` | 1.0 | 0.9999999798763621 | 0.9999999787878232 | 2.0e-8 |
| 79372.2 | 0.106553 | 0.9999999999999997 | `cdf` | 1.0 | 0.921339739174892 | 0.9213397391748921 | 8.5e-2 |
| 1e4 | 1e4 | 0.5 | `cdf` | 0.5000000000286743 | 0.5 | 0.4999999999999989 | 5.7e-11 |
| 1e4 | 1e4 | 0.46464554479100484 | `cdf` | 6.739002403522573e-24 | 6.739002403910842e-24 | 6.739002403910101e-24 | 5.8e-11 |
| 1 | 1e4 | 1e-100 | `cdf` | 9.99999999975849e-97 | 1.0000000000000001e-96 | 1.0000000000000001e-96 | 2.4e-11 |
| 1e4 | 0.1 | 0.9998951379361433 | `cdf` | 0.02233056971103642 | 0.022330569703958216 | 0.02233056970395821 | 3.2e-10 |
| 0.1 | 1e4 | 0.00010486206385669913 | `sf` | 0.022330569715089843 | 0.022330569703971435 | 0.022330569703971435 | 5.0e-10 |
| 71839.6 | 0.131114 | 0.9999877675 | `cdf` | 0.03949630795747194 | 0.039496308017432535 | 0.03949630801743255 | 1.5e-9 |

Second opinion: SciPy (`scipy.stats.beta`) agrees with mpmath within
2.2e-15 at every row except two: 1.1e-13 at a = b = 1e4,
x = 0.46464554479100484 (within tolerance), and a = b = 0.5,
x = 1 - 1e-15, where SciPy is 1.1e-9 off itself. There the closed form I = (2/pi) asin(sqrt x) gives the
mpmath value (tested in `reference/tests/test_functions.py`, together with
the positive-term series DLMF 8.17.8 at the 3rd, 5th, 7th and 8th rows).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double (one `beta_cdf` point, a = 1, b = 100,
   x = 1e-100, is an exact rounding tie and tagged `unstable-reference`).
   The route, the continued fraction DLMF 8.17.22 on whichever side of the
   switch point it converges with 1 - x exact (`_mp_beta_cdf`,
   `_mp_beta_sf`), is tested against `mpmath.betainc`, the positive-term
   series DLMF 8.17.8 in both tails, closed forms for a = b = 1/2, b = 1 and
   a = 1, and cdf + sf = 1, in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random points (Python `random.Random(20261015)`):
   a, b = 10^u with u uniform in [-1, 5]; x uniform in (0, 1) for half the
   points, x = 10^v with v uniform in [-300, 0] for a quarter, and
   x = 1 - 10^v with v uniform in [-16, 0] for the rest; each a, b and 10^v
   rounded to 6 significant digits; mpmath at 50 and 100 digits (they agreed
   at all 20,000). `cdf` fails 929 (4.6%), `sf` 910 (4.6%). Removing the
   cut-off (cause 1) fixes 45 `cdf` points, all with 1 - x <= 1e-15 and
   b < 0.98; the worst `cdf` error drops from 8.5e-2 to 1.5e-9. Every other
   failure has max(a, b) >= 10.3; worst per decade of max(a, b), `cdf`:
   2.7e-12 (10 to 100), 4.8e-12 (100 to 1000), 4.3e-11 (1000 to 1e4),
   1.5e-9 (1e4 to 1e5); `sf`: 1.6e-12, 2.4e-12, 4.0e-11, 2.7e-10.
5. **Root cause** (confirmed with a copy of `beta_reg` and of `Beta::cdf`
   and `Beta::sf`, scratch, not committed, that reproduces statrs bit for
   bit on all 1,732 interior grid points and the 20,000 random points and
   can switch the pieces below). `beta_reg` (`src/function/beta.rs`) is the
   Numerical Recipes continued fraction (at most 140 iterations) with the
   prefactor `bt = exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) +
   a ln x + b ln(1 - x))`; for x >= (a + 1)/(a + b + 2) it returns
   `1 - bt * cf(b, a, 1 - x) / b`. `Beta::sf` computes
   `1 - beta_reg(a, b, x)` below that switch point and
   `beta_reg(b, a, 1 - x)` above it.
   - **Cause 1:** `beta_reg` sets `bt = 0` when x is within 1e-15 of 1
     (`ulps_eq!(x, 1.0, epsilon = MODULE_EPS)`, `MODULE_EPS = 1e-15`), so it
     returns 1 - 0 = 1. The neglected tail is
     I_(1-x)(b, a), about (1 - x)^b / (b B(a, b)), which for b < 1 is far
     from negligible: (1e-15)^0.1 = 0.03. `sf` is not affected (it calls
     `beta_reg(b, a, 1 - x)` there, where 1 - x is not near 1). Removing the
     cut-off removes all 19 grid points of group 1.
   - **Cause 2:** the prefactor's logarithms are rounded separately and
     then summed and exponentiated. They reach 1.78e5 on the grid
     (`ln_gamma(2e4)` at a = b = 1e4), so rounding alone can cost
     2.0e-11 (`reference/tests/test_functions.py` checks this number), and
     statrs's `ln_gamma` has its own error of about 1e-13 per call
     (finding 010). Using the correctly rounded prefactor (computed with
     mpmath) instead, with the cut-off also removed, leaves **no failures
     on the grid** (worst then 7.2e-13 for `cdf`, 7.9e-13 for `sf`) and
     8 of 20,000 off it (7 `cdf`, 1 `sf`, 1.1e-12 to 7.4e-12, all with
     a or b >= 15,576 and x within 1e-3 of 0 or 1: the continued
     fraction's own rounding).
   - **Amplification:** where the result is 1 minus the other tail
     (`beta_reg` above the switch point, `sf` below it) and that tail is
     not small, its relative error is multiplied by the ratio of the two
     tails: at a = 0.1, b = 1e4, x = 0.00010486 `cdf` is 0.978 with an
     error of 1.14e-11, and `sf = 1 - cdf` = 0.0223 has 44 times that,
     5.0e-10 (`cdf` at the mirrored point a = 1e4, b = 0.1 likewise). Moving the switch
     to the mean a / (a + b), as finding 001's fix does, is no remedy here:
     the continued fraction then converges too slowly within 140 iterations
     (grid failures rise to 10 even with the exact prefactor; worst 9.1e-7).
   - The 140-iteration cap alone changes nothing on the grid.

## Notes

Boost.Math's `ibeta`/`ibetac` compute the prefactor with a Lanczos
approximation that combines the gamma terms before rounding, have no
cut-off near 1, and return whichever of I and 1 - I is small directly.
