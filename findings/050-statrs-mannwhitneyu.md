# 050: statrs mannwhitneyu returns the wrong tail for exact p-values when x is not the longer sample, and small asymptotic p-values as 1 - cdf

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::mannwhitneyu::mannwhitneyu`)
- **Severity:** high (silently wrong answers: SciPy's documented example,
  females against males with the exact method and "less", gives p = 0.944
  instead of 0.0556; three values all below six others give a two-sided p
  of 1 instead of 0.0238; a two-sided asymptotic p-value of 7.1e-18 is
  returned as 0)
- **Status:** approved (Claude, checklist complete, 2026-10-29)
- **Affected inputs:** 180 of 690 points of `reference/data/mannwhitneyu.jsonl`
  (tolerance 1e-13, for the statistic and the p-value). The exact list is
  `MANNWHITNEYU` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes
  (section "Root cause"):
  1. **Exact p-value complemented when n1 <= n2** (120 points: every exact
     point, method `Exact` or `Automatic`, whose x is not longer than y; 4.1%
     to 3.1e6 times off). statrs returns 1 - P(U >= u), i.e. P(U < u),
     instead of P(U >= u): small p-values come out near 1 and large ones near
     0, and two-sided ones are clipped to 1. With x longer than y (all 90 such
     exact points) the answer is right. The statistic U1 is always right.
  2. **The asymptotic p-value as 1 - cdf** (60 points of methods
     `Automatic`, `AsymptoticInclContinuityCorrection` and
     `AsymptoticExclContinuityCorrection`, up to 1.0e-10 off on this grid,
     whose p-values are at least 3.2e-7): `1.0 - norm_dist.cdf(z)` loses
     digits as p falls (6.8e-7 off at p = 1.4e-11) and is 0 below about
     1e-16, as in finding 026.

  Another 268 points of the dataset fail only through statrs's normal
  cdf, whose `erfc` is off by up to 1.03e-10 (finding 009); they are listed
  under finding 009 as `MANNWHITNEYU_ERFC`.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::mannwhitneyu::{MannWhitneyUMethod as M, mannwhitneyu};
use statrs::stats_tests::Alternative as A;

fn main() {
    let r = |x: &[f64], y: &[f64], m: M, a: A| mannwhitneyu(x, y, m, a).unwrap();
    let range = |a: u32, b: u32| (a..=b).map(f64::from).collect::<Vec<_>>();
    // SciPy's documented example: females vs males, exact, "less": p = 0.0556
    let males = [19.0, 22.0, 16.0, 29.0, 24.0];
    let females = [20.0, 11.0, 17.0, 12.0];
    println!("{:?}", r(&females, &males, M::Exact, A::Less));
    println!("{:?}", r(&males, &females, M::Exact, A::Greater)); // the same test, x longer: right
    // Every x below every y: two-sided p = 2 / C(9, 3) = 0.0238
    println!("{:?}", r(&range(1, 3), &range(4, 9), M::Exact, A::TwoSided));
    // Automatic picks the exact method (n1 = 8): "less" p = 1 / C(30, 8) = 1.7e-7
    println!("{:?}", r(&range(1, 8), &range(9, 30), M::Automatic, A::Less));
    // Normal approximation, 50 vs 50 values, fully separated
    println!("{:?}", r(&range(1, 50), &range(51, 100), M::AsymptoticInclContinuityCorrection, A::TwoSided));
    println!("{:?}", r(&range(51, 80), &range(1, 30), M::AsymptoticExclContinuityCorrection, A::Greater));
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below.

## Expected

| x, y, method, alternative | U1 (all agree) | p, crate | p, mpmath (true value) | p, SciPy |
|---|---|---|---|---|
| females, males, exact, less | 3 | 0.9444444444444444 | 0.05555555555555555 | 0.05555555555555555 |
| males, females, exact, greater | 17 | 0.05555555555555555 | 0.05555555555555555 | 0.05555555555555555 |
| 1..3, 4..9, exact, two-sided | 0 | 1.0 | 0.023809523809523808 | 0.023809523809523808 |
| 1..8, 9..30, automatic, less | 0 | 0.9999998291452564 | 1.7085474356838674e-7 | 1.7085474356838674e-7 |
| 1..50, 51..100, asymptotic with correction, two-sided | 0 | 0.0 | 7.066071930388956e-18 | 7.066071930388932e-18 |
| 51..80, 1..30, asymptotic without correction, greater | 900 | 1.4359735622804237e-11 | 1.4359745331601645e-11 | 1.4359745331601617e-11 |

The second row is the same test with the samples swapped, which statrs gets
right. Second opinion: SciPy (`scipy.stats.mannwhitneyu`, which computes
the exact distribution by a different recurrence) agrees with mpmath within
2.1e-15 at every point of the dataset and every row above. The exact rows
need no numerics: 0.0556 = 7/126 is SciPy's documented value (7 of the
C(9, 4) = 126 splits have U2 >= 17), and with every x below every y only
one split of C(9, 3) = 84 and one of C(30, 8) = 5852925 is as extreme.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route (`_mp_mannwhitneyu`)
   takes U1 from exact average ranks (rational arithmetic), the exact
   p-value from the null counts of the recurrence c(m, n; u) = c(m - 1, n;
   u - n) + c(m, n - 1; u) (Mann and Whitney 1947) as an exact fraction,
   and the asymptotic one from erfc with the tie-corrected variance. New
   tests in `reference/tests/test_functions.py` check the counts against
   enumeration of every split, U1 against pair counting, the tie-corrected
   variance against the exact variance over every split of tied samples,
   the erfc tail against `mpmath.ncdf`, the automatic rule, and SciPy's
   published example (7/126).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261029)`):
   n1 and n2 uniform in 1..30; integer values uniform in 0..5, 0..50 or
   0..10^6 (the small ranges give ties), x shifted by a random integer, 30%
   of samples rescaled by 1e-3, 7.25 or 1e5; method and alternative uniform;
   exact points only without ties and with at most 5 million splits; mpmath
   at 50 and 100 digits (they agreed at all points); points whose rounding
   bound (the tolerance reason's) exceeds the tolerance left out. Of the
   220 exact points with n1 <= n2, 216 fail (4.4% to 1.6e6 times off); the
   4 that pass have a true p of 1 or 0.5, where the complement happens to
   agree. None of the 187 exact points with n1 > n2 fails. 328 asymptotic
   points fail through cause 2 (1.8e-11 to 1.7e-10 off) and 1,454 more only
   through finding 009's `erfc` (up to 1.03e-10).
5. **Root cause** (from the source, `src/stats_tests/mannwhitneyu.rs`,
   confirmed with a copy that is bit for bit equal to statrs on all 690 grid
   points and all 3,000 random ones):
   - **Cause 1:** `calc_mwu_exact_pvalue` enumerates every split of the
     ranks into a group of k = min(n1, n2) and counts those whose U is at
     least u, which is P(U >= u) whichever sample has k values (the null
     distribution of U is the same for both). It then returns
     `1.0 - numerator as f64 / total as f64` when `k == n1` and the plain
     ratio otherwise; the comment says the complement is for the other
     branch. Dropping the `1.0 -` makes every exact grid and random point
     pass.
   - **Cause 2:** `1.0 - norm_dist.cdf(z)` in
     `calc_mwu_asymptotic_pvalue`. Using the upper tail directly
     (`0.5 * erfc(z / sqrt 2)`) removes this cause; the p-values that then
     still fail are exactly the `erfc` error of finding 009: the copy with
     both remedies and an accurate `erfc` (libm 0.2.16) passes every grid
     and random point.

## Notes

Cause 1 hits the common case: with the exact method (chosen automatically
when a sample has at most 8 values and there are no ties), any call whose
first sample is the shorter one, or as long as the other, gets the wrong
tail. A one-sided test then reports p near 1 for a clear effect, or a tiny
p where there is none.

Also seen, not on the grid (no listed points): the exact method enumerates
all C(n1 + n2, min(n1, n2)) splits, and `Automatic` chooses it whenever one
sample has at most 8 values. Measured in release mode: 0.33 s for 8 and 30
values, 2.5 s for 8 and 40; at the same rate 8 and 100 values would take
about 40 minutes (an estimate, not run). SciPy's recurrence needs
milliseconds. The grid keeps the exact sizes at most 3.2 million splits.
