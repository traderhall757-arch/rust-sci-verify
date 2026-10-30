# 052: statrs Hypergeometric pmf is NaN for populations above about 1030, and its cdf loses digits

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Hypergeometric`:
  `pmf` and `cdf`, seen through `statrs::stats_tests::fisher::fishers_exact`)
- **Severity:** high (silently wrong answers through it: a two-sided
  Fisher p-value of 0.493 instead of 0.051; one-sided p-values up to
  6.1e-12 off. `pmf` itself is NaN at the most probable value for a
  population of 1900)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 114 of 411 points of `reference/data/fisher_exact.jsonl`
  (tolerance 1e-13), which pass with an accurate pmf and cdf and statrs's
  Fisher code otherwise unchanged.
  The exact list is `FISHER_EXACT_HYPERGEOMETRIC` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). `Hypergeometric`
  has no dataset of its own yet. Two causes (section "Root cause"):
  1. **pmf overflow** (8 points, all two-sided; 5 more grid points have it
     together with finding 051's cause 1 and are listed there): `pmf` is a
     ratio of binomial coefficients in double precision, and a coefficient or
     the product of the top two is infinite once it exceeds 1.8e308, so the
     pmf is NaN (inf / inf, or inf x 0 outside the support) or 0 for
     populations from about 1030 values. The two-sided p-value compares
     table probabilities and then returns an arbitrary number: 0.493
     instead of 0.051, 4.8e-13 instead of a value below the smallest double.
  2. **cdf terms from logarithms** (106 points: all 80 one-sided failures, 1e-13
     to 1.5e-12 off on the grid, and 26 two-sided ones; 14 of these also
     pass with finding 051's remedy instead): each term is
     exp(ln C(K, x) + ln C(N - K, n - x) - ln C(N, n)), and the rounding
     error of logarithms of size up to ln(N!) (about 7,300 to 12,400 here) becomes a relative error of the
     term. One-sided failures start at a total of 182 on the grid (135 off
     it).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Discrete, DiscreteCDF, Hypergeometric};
use statrs::stats_tests::fisher::fishers_exact;
use statrs::stats_tests::Alternative as A;

fn main() {
    let f = |t: [u64; 4], a: A| fishers_exact(&t, a).unwrap();
    let h = Hypergeometric::new(1900, 1000, 700).unwrap();
    println!("{:?} {:?}", h.pmf(0), h.pmf(368));
    println!("{:?}", f([0, 1000, 700, 200], A::TwoSided));
    println!("{:?}", f([389, 611, 311, 589], A::TwoSided));
    println!("{:?}", Hypergeometric::new(1200, 600, 300).unwrap().cdf(90));
    println!("{:?}", f([90, 510, 210, 390], A::Less));
    println!("{:?}", f([173, 1902, 317, 193], A::Less));
}
```

These are the last lines of the program that was run (the first ones are
in finding 051; release and debug mode give the same output); it gives the
"Crate" column below.

## Expected

| Call | Cause | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|
| pmf(0), N = 1900, K = 1000, n = 700 | - | 0.0 | 0.0 (about 1.7e-336) | 0.0 |
| pmf(368), same | 1 | NaN | 0.03795160236675016 | 0.03795160236675012 |
| Fisher [0, 1000, 700, 200], two-sided | 1 | 4.759526106568046e-13 | 0.0 (below 5e-324) | 0.0 |
| Fisher [389, 611, 311, 589], two-sided | 1 | 0.49280766573914747 | 0.05100887098887347 | 0.05100887098887342 |
| cdf(90), N = 1200, K = 600, n = 300 | 2 | 5.382263667531292e-16 | 5.382263667523354e-16 | 5.382263667523352e-16 |
| Fisher [90, 510, 210, 390], less | 2 | 5.382263667531292e-16 | 5.382263667523354e-16 | 5.382263667523352e-16 |
| Fisher [173, 1902, 317, 193], less (off the grid) | 2 | 1.3668570050659537e-141 | 1.366857005074232e-141 | 1.3668570050742307e-141 |

The first row is right (the true value underflows) and shows that the
NaN is not an underflow. Second opinion: SciPy (`scipy.stats.hypergeom`
and `fisher_exact`) agrees with mpmath within 1.1e-15 at every row above
and within 2.1e-15 at every point of the dataset. The pmf and cdf rows are
exact fractions of binomial coefficients, computed in integer arithmetic.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: as in finding 051 (exact table counts, tested
   against enumeration of labellings and against `mpmath.binomial` at 60
   digits for the largest margins).
3. **Second opinion**: see "Expected".
4. **Off-grid**: the 3,000 random tables of finding 051 (seed 20261030).
   954 fail only through `Hypergeometric`: all 714 failing one-sided points
   (cause 2, 1e-13 to 6.1e-12 off, totals from 135) and 240 two-sided ones,
   38 of which need the pmf remedy (relative errors up to 7.1e217; 3 of
   these have a finite C(N, n), e.g. (1400, 32, 190, 10): statrs's
   search starts at 0, outside the support (1390..1432), where C(1432, x)
   is infinite and C(200, 1590 - x) is 0, so the pmf is inf x 0 = NaN).
5. **Root cause** (from the source, `src/distribution/hypergeometric.rs` and
   `src/function/factorial.rs`, confirmed with the copy of finding 051, bit
   for bit equal to statrs on all grid and random points):
   - **Cause 1:** `pmf` returns `binomial(K, x) * binomial(N - K, n - x) /
     binomial(N, n)`, and `factorial::binomial` is
     `floor(0.5 + exp(ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)))`,
     +inf once the coefficient exceeds the largest double (C(1030, 515) is
     about 2.9e308). Taking the pmf as the exponential of the log terms that
     `cdf` uses removes the NaN; of the 8 points, 2 then pass and 6 also
     need the cause 2 remedy (5 more need finding 051's remedy as well).
   - **Cause 2:** `cdf` sums exp(ln_binomial(K, i) + ln_binomial(N - K, n - i)
     - ln_binomial(N, n)) for i = 0..=x. The logarithms are large (ln C(1900,
     700) = 1246) and `ln_factorial` above 170 comes from `ln_gamma`, so
     each exponent carries an absolute rounding error of order 1e-12, which
     exp turns into the same relative error. Computing the terms instead by
     the ratio P(x + 1) / P(x) = (K - x)(n - x) / ((x + 1)(N - K - n + x + 1))
     outward from the mode, normalised by their sum, makes every one-sided
     grid point pass, with statrs's pmf left as it is.

   With both remedies and finding 051's, the copy passes every grid and
   random point (worst 3.5e-15).

## Notes

`Hypergeometric::sf` and `ln_pmf` use the same log terms as `cdf` (not
checked; no listed points). A dataset for the hypergeometric distribution
itself (with the other discrete distributions) is not on the roadmap yet;
the owner decides whether it belongs in Milestone 4.
