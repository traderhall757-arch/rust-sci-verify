# 051: statrs fishers_exact gives wrong two-sided p-values: an upper tail taken as 1 - cdf, and a 1e-4 tie window

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::fisher::fishers_exact`,
  and `fishers_exact_with_odds_ratio`, which calls it)
- **Severity:** high (silently wrong answers: a negative p-value,
  -2.1e-13, for a table of statrs's own tests whose true p is 9.0e-16;
  2.2e-14 instead of 2.2e-59; 0.529 instead of 0.368 for the 6 x 9 x 16 x
  41 table)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 51 of 414 points of `reference/data/fisher_exact.jsonl`
  (tolerance 1e-13, for the odds ratio and the p-value), all two-sided.
  The exact list is `FISHER_EXACT` in `crates/check-statrs/src/expected.rs`,
  the report is [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md).
  The odds ratio is always right, and so are the one-sided p-values apart
  from finding 052. Two causes (section "Root cause"):
  1. **Upper tail as 1 - cdf** (42 points, 2.9e-12 to 8.1e273 times off):
     the two-sided p-value adds an upper tail computed as `1 - cdf`, which
     keeps only an absolute accuracy of about 1e-16 times the number of
     terms. Small p-values lose digits (3e-4 relative error at 1.45e-11), and
     below about 1e-14 the result is rounding noise: 2.2e-14 where the true
     value is 2.2e-59, and even negative.
  2. **Tie window of 1e-4** (12 points, 38% to 91% off; 3 of them also have
     cause 1): a table counts as "at most as probable as the observed one"
     if its probability is up to 1.0001 times the observed one, so tables
     that are more probable are added to the p-value. SciPy uses a window of
     1e-14, which only absorbs rounding.

  Another 114 points of the dataset fail only through statrs's
  `Hypergeometric` pmf and cdf (finding 052); they are listed there as
  `FISHER_EXACT_HYPERGEOMETRIC`.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::fisher::fishers_exact;
use statrs::stats_tests::Alternative as A;

fn main() {
    let f = |t: [u64; 4], a: A| fishers_exact(&t, a).unwrap();
    println!("{:?}", f([172, 46, 90, 127], A::TwoSided)); // statrs's own test table
    println!("{:?}", f([0, 20, 20, 0], A::TwoSided));
    println!("{:?}", f([0, 100, 100, 0], A::TwoSided));
    println!("{:?}", f([6, 9, 16, 41], A::TwoSided));
    println!("{:?}", f([1, 28, 32, 25], A::TwoSided));
}
```

These are the first lines of the program that was run (release and debug
mode give the same output; the rest of it is in finding 052); it gives the
"Crate" column below.

## Expected

| Table [a, b, c, d], two-sided | Cause | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|
| 172, 46, 90, 127 | 1 | -2.1425722797266258e-13 | 9.041009036526075e-16 | 9.041009036526073e-16 |
| 0, 20, 20, 0 | 1 | 1.4504508705215358e-11 | 1.4508889103849688e-11 | 1.4508889103849688e-11 |
| 0, 100, 100, 0 | 1 | 2.2315482794965646e-14 | 2.2087606931995028e-59 | 2.208760693199503e-59 |
| 6, 9, 16, 41 | 2 | 0.5294033745950272 | 0.3676007792428934 | 0.36760077924289336 |
| 1, 28, 32, 25 | 2 | 9.27084059632044e-7 | 4.858505327102482e-7 | 4.858505327102481e-7 |

The true values are exact fractions: [0, 20, 20, 0] has p = 2 / C(40, 20)
(only it and [20, 0, 0, 20] are as improbable), [0, 100, 100, 0] has
p = 2 / C(200, 100). Second opinion: SciPy (`scipy.stats.fisher_exact`)
agrees with mpmath within 2.1e-15 at every point of the dataset and within
2.2e-16 at every row above. statrs's own test passes on the first table
because `prec::assert_relative_eq!` also accepts any absolute difference
below `DEFAULT_EPS` = 1e-9.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The route
   (`_mp_fisher_exact`) counts tables exactly: C(n1, x) C(n2, n - x) in
   integer arithmetic, compared exactly for the two-sided sum, divided by
   C(n1 + n2, n) once. New tests in `reference/tests/test_functions.py`
   check the counts against enumeration of every labelling of small
   tables, "greater" against "less" on the second column, all three
   alternatives against `mpmath.binomial` probabilities at 60 digits for
   the largest margins, the published values (SciPy's 5/143 and the
   (6, 2, 1, 4) example, Fisher's tea tasting 17/70 and 34/70), and that
   each near-tie table has a more probable table within 1e-4 that the
   p-value leaves out.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random tables (Python `random.Random(20261030)`):
   total log-uniform in 4..3000, n1 and n uniform, a uniform in its support
   or at the mean + k sd with k uniform in (-10, 10), b c > 0, alternative
   uniform; true p from the exact counts. 370 of the 1,032 two-sided points
   fail through cause 1 (all still fail with finding 052 remedied). No
   random table had a near tie; a search over all margins (totals 72 to
   117, the first 40 tables with another table at most 1e-4 more probable)
   gave 40 failures through cause 2 (38% to 91% off) that remain with
   every other remedy.
5. **Root cause** (from the source, `src/stats_tests/fisher.rs`, confirmed
   with a copy that is bit for bit equal to statrs on all 414 grid points
   and all 3,000 random ones):
   - **Cause 1:** `1.0 - dist.cdf(table[0] - 1)` (upper tail when a is at
     or above the mode) and `p_lower + 1.0 - dist.cdf(guess - 1)` (the other
     tail when a is below it). Summing the upper tail directly, as
     `Hypergeometric::sf` does, removes this cause.
   - **Cause 2:** `const EPSILON: f64 = 1.0 - 1e-4`, used both for the
     "p = 1 at the mode" test and in `binary_search`, which accepts
     `pmf(guess) <= p_exact / EPSILON`. With 1 - 1e-14 (SciPy's window) all
     12 grid points and the 40 searched tables pass.

   The copy with both remedies and finding 052's passes every grid and
   random point (worst 3.5e-15).

## Notes

Cause 1 hits the case a two-sided test is used for: a clear effect. The
smaller the two-sided p-value, the fewer digits are right; below about
1e-14 none are, and the result can be negative.

Not on the grid (no listed points): tables with an empty row or column
return p = 1 and odds ratio NaN by convention, as in SciPy, and
`fishers_exact_with_odds_ratio` forms a d and b c in `u64`, which overflows
for products above 1.8e19 (a panic in debug mode; not run). Tables with b c = 0 (odds ratio +inf) are left out
of the grid because the generator would tag them `overflow`, which leaves
the p-value unjudged.
