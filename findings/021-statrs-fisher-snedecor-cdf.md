# 021: statrs FisherSnedecor cdf and sf round their beta argument: sf is 0 in the upper tail, both return exactly 1 and lose digits

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::FisherSnedecor`,
  methods `cdf` and `sf` of `ContinuousCDF`; both call
  `statrs::function::beta::beta_reg`)
- **Severity:** high (silently wrong answers: `sf` returns 0 for every x
  above about 1e16 d2 / d1, e.g. 0 for 1e-20 at d1 = 5, d2 = 2, x = 1e20;
  `cdf` returns exactly 1 for 0.99753 at d1 = 1e5, d2 = 0.5, x = 1e10;
  `sf` is 2.6% off at that point and up to 138% off the grid), plus
  finding 019's accuracy loss for large degrees of freedom (up to 1.3e-9)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 194 of 1999 points of `reference/data/f_cdf.jsonl`
  and 440 of 1999 points of `reference/data/f_sf.jsonl` (tolerance 1e-12);
  the exact lists are `F_CDF` and `F_SF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). By symptom
  (causes in "Root cause"):
  1. **`sf` exactly 0 in the upper tail** (117 `sf` points: every x = 1e20,
     1e100, 1e300 where the true value is a normal double, d1 x / d2 at
     least 1.7e18): the true value is up to 7.8e-6 (d1 = 1e5, d2 = 0.5,
     x = 1e20). Any x with d1 x / d2 above about 2^53 = 9.0e15 gives 0.
  2. **Exactly 1 next to the ends** (20 `cdf` and 20 `sf` points, all with
     the other degree of freedom 0.5 or 1): `cdf` at x = 1e20 (and 1e10 for
     d1 = 1e5) and its mirror image `sf` at x = 1e-20 (and 1e-10 for
     d2 = 1e5) return 1 where the true value is 1 - 5.4e-11 up to
     1 - 2.5e-3 (d1 = 1e5, d2 = 0.5, `cdf(1e10)` = 0.99753).
  3. **Digits lost from the rounded argument** (the rest of the points
     with max(d1, d2) < 100: 8 `cdf`, 72 `sf`, plus some with large
     degrees of freedom): from 1e-12 up to 1.3e-6 on the grid
     (`sf(1e10)` at d1 = 5, d2 = 30), growing with d1 x / d2 for `sf` and
     `cdf` above the median, and with d2 / (d1 x) for `sf` below it.
     With large d1 they reach 2.6e-2 (`sf(1e10)` = 0.002532 instead of
     0.002467 at d1 = 1e5, d2 = 0.5).
  4. **Large degrees of freedom** (the remaining points, every one with
     max(d1, d2) >= 100; 150 `cdf` points, most of the remaining `sf`
     points): errors from 1e-12 to 1.3e-9 across the bulk and both tails
     (`cdf(1)` = 0.4999999993270384 at d1 = d2 = 1e5): finding 019's
     prefactor, plus the 140-iteration cap of the continued fraction.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, FisherSnedecor};

fn main() {
    let f = |d1, d2| FisherSnedecor::new(d1, d2).unwrap();
    println!("{:?}", f(5.0, 2.0).sf(1e20)); // should be 1e-20
    println!("{:?}", f(1e5, 0.5).sf(1e10)); // should be 0.0024669656598050476
    println!("{:?}", f(1000.0, 30.0).sf(1e10)); // should be 4.122892850626986e-145
    println!("{:?}", f(10.0, 0.5).cdf(1e20)); // should be 0.99999234586564
    println!("{:?}", f(1e5, 0.5).cdf(1e10)); // should be 0.9975330343401949
    println!("{:?}", f(0.5, 1e5).sf(1e-10)); // should be 0.9975330343401949
    println!("{:?}", f(1000.0, 0.5).cdf(1e10)); // should be 0.9975334923052093
    println!("{:?}", f(1e5, 1e5).cdf(1.0)); // should be 0.5
    println!("{:?}", f(1e5, 1.0).cdf(0.5)); // should be 0.1573023203491649
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| d1 | d2 | x | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|---|
| 5 | 2 | 1e20 | `sf` | 0.0 | 1e-20 | 1e-20 | 1 |
| 1e5 | 0.5 | 1e10 | `sf` | 0.002532302953833911 | 0.0024669656598050476 | 0.002466965659805048 | 2.6e-2 |
| 1000 | 30 | 1e10 | `sf` | 4.123813529802137e-145 | 4.122892850626986e-145 | 4.122892850626836e-145 | 2.2e-4 |
| 10 | 0.5 | 1e20 | `cdf` | 1.0 | 0.99999234586564 | 0.99999234586564 | 7.7e-6 |
| 1e5 | 0.5 | 1e10 | `cdf` | 1.0 | 0.9975330343401949 | 0.9975330343401949 | 2.5e-3 |
| 0.5 | 1e5 | 1e-10 | `sf` | 1.0 | 0.9975330343401949 | 0.9975330343401949 | 2.5e-3 |
| 1000 | 0.5 | 1e10 | `cdf` | 0.9975339853092353 | 0.9975334923052093 | 0.9975334923052094 | 4.9e-7 |
| 1e5 | 1e5 | 1 | `cdf` | 0.4999999993270384 | 0.5 | 0.49999999999999867 | 1.3e-9 |
| 1e5 | 1 | 0.5 | `cdf` | 0.15730232041574932 | 0.1573023203491649 | 0.157302320349165 | 4.2e-10 |

Second opinion: SciPy (`scipy.stats.f`) agrees with mpmath within 2.7e-15
at every row except the third, where it is 3.6e-14 off (within
tolerance); its `sf` forms d2 / (d2 + d1 x) directly. Where SciPy is
itself wrong on the grid (6 `cdf`, 11 `sf` points, up to 6.4e-11), the
closed form F = (2/pi) atan(sqrt x) for d1 = d2 = 1 and the positive-term
series DLMF 8.17.8 confirm mpmath (`reference/tests/test_functions.py`,
which also checks rows 2, 5, 6, 8 against the series).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double (they agree at every point). The route
   (`_mp_f_cdf`, `_mp_f_sf`: the continued fraction DLMF 8.17.22 on whichever
   side of the switch point it converges, with z = d1 x / (d1 x + d2) and
   y = d2 / (d1 x + d2) formed separately at the working precision) is tested
   against `mpmath.betainc`, the positive-term series on the smaller tail,
   closed forms for d1 = d2 = 1, d2 = 2 and d1 = 2, and cdf + sf = 1, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random points (Python `random.Random(20261015)`):
   d1, d2 = 10^u with u uniform in [-0.5, 5]; for half the points
   x = exp(2 s g) with g standard normal and s = sqrt(2/d1 + 2/d2) (the
   bulk), for the rest x = 10^v with v uniform in [-300, 300]; d1, d2 and x
   rounded to 6 significant digits; mpmath at 50 and 100 digits (they agreed
   at all 20,000). `cdf` fails 3,601 (18%), worst 1.4e-3 (exactly 1 at
   d1 = 50.1219, d2 = 0.405702, x = 3.1148e13); `sf` fails 4,421 (22%),
   worst 1.38 (2.4 times the true value at d1 = 1075.65, d2 = 47.2902,
   x = 2.16193e13), and returns 0 for a positive true value at 1,019 of them.
   Every failure not explained by causes 1 and 2 has max(d1, d2) >= 22.8.
   With all four causes below switched off (copy described in item 5), 26
   `cdf` and 31 `sf` failures remain, 1.0e-12 to 5.5e-12, all with one
   degree of freedom above 3.1e4 and the other at most 185: the continued fraction's own rounding,
   amplified where the result is 1 minus the larger tail (finding 019,
   "Amplification").
5. **Root cause** (confirmed with a copy of `beta_reg` and of
   `FisherSnedecor::cdf` and `sf`, scratch, not committed, that reproduces
   statrs bit for bit on all 1999 grid points and the 20,000 random points
   and can switch the pieces below). `cdf` is
   `beta_reg(d1/2, d2/2, z)` and `sf` is `beta_reg(d2/2, d1/2, 1 - z)`, with
   `z = d1 x / (d1 x + d2)` (`src/distribution/fisher_snedecor.rs`).
   - **Cause 1, the argument:** 1 - z is formed from the rounded z. Its
     relative error is about 2^-53 z / (1 - z) = 2^-53 d1 x / d2, unbounded
     as x grows; once d1 x / d2 exceeds about 2^53, z rounds to 1 and
     `sf` = `beta_reg(b, a, 0)` = 0 (symptom 1). Inside `beta_reg` the same
     happens for `cdf` above the switch point (it computes 1 - z again), and
     for `sf` at small x, where 1 - z is close to 1 and `beta_reg` uses
     1 - (1 - z). Passing y = d2 / (d1 x + d2), formed directly, instead
     of 1 - z (and to `beta_reg` instead of 1 - x) removes 22 `cdf` and 241
     `sf` grid failures, including every exact 0.
   - **Cause 2, finding 019's cut-off:** `beta_reg` sets its prefactor to 0
     when its argument is within 1e-15 of 1 (`MODULE_EPS`), so it returns
     exactly 1 (symptom 2): `cdf` once d2 / (d1 x) < 1e-15, `sf` once
     d1 x / d2 < 1e-15. The neglected tail is about C (1e-15)^(d/2), far
     from negligible for d below 2. Removing the cut-off as well removes 20
     more points of each.
   - **Cause 3, finding 019's prefactor:** the logarithms in
     `exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a ln z + b ln y)`
     are rounded separately; they reach 1.05e6 here (ln Gamma(1e5)), so the
     rounding alone can cost 1.2e-10 (tested). With the correctly rounded
     prefactor and arguments (computed with mpmath) as well, 7 `cdf` and
     13 `sf` grid points remain.
   - **Cause 4, the iteration cap:** the continued fraction stops after 140
     iterations, too few near the median for a, b around 5e4:
     `cdf(1)` at d1 = d2 = 1e5 keeps an error of 1.2e-9 with
     causes 1 to 3 removed, and 0 without the cap (finding 019's grid had
     shapes only up to 1e4, where the cap changed nothing). Lifting it
     leaves 6 `cdf` and 12 `sf` points, 1.1e-12 to 3.6e-12, all with one
     degree of freedom 1e5 and the other at most 30, results from 1e-4 to
     0.15 formed as 1 minus the larger tail (the residue described in
     item 4).

## Notes

SciPy's `fdtr`/`fdtrc` (Cephes, then Boost's `ibeta`) form
w = d1 x / (d2 + d1 x) for the cdf and w = d2 / (d2 + d1 x) for the
survival function, so neither tail rounds 1 - z; SciPy's remaining errors
on this grid (cdf at most 1.7e-12 with d1 = 1e5) come from the cdf's z near
1. `FisherSnedecor::inverse_cdf` has its own problems:
[finding 022](022-statrs-fisher-snedecor-inverse-cdf.md).
