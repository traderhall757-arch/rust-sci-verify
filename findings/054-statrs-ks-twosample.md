# 054: statrs ks_twosample never returns for equal empirical cdfs, and gives wrong statistics at D = 0 and wrong exact and asymptotic p-values

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::ks_test::ks_twosample`,
  all four `KSTwoSampleAlternativeMethod`s)
- **Severity:** high (the call never returns when the two samples have the
  same empirical cdf; a negative two-sided exact p-value, -6.0e-15 instead
  of 2.2e-59)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 30 of 363 points of `reference/data/ks_2samp.jsonl`
  (tolerance 5e-12, for the statistic and the p-value). The exact list is
  `KS_2SAMP` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Hodges's
  one-sided approximation itself is computed right. Four causes (section
  "Root cause"):
  1. **Statistic not 0 where it should be** (12 points, methods 0 and 1):
     the empirical cdfs are built by adding 1/n1 and 1/n2 once per value,
     so they end slightly off 1 (for example 0.9999999999999999) and their
     difference is rounding noise where it should be exactly 0. A one-sided
     statistic whose true value is 0 (D- when all of x is below y, for
     example) comes out as 5.6e-17 to 3.3e-16, and its p-value as
     0.9999999999999997 instead of 1. With ties, the two-sided D = 0 is
     affected the same way. In large samples (sizes near 1e5, off the grid)
     the long sums also move statistics above 0, by up to 4.0e-11.
  2. **Never returns** (4 points, method 3, samples with equal empirical
     cdfs, D = 0): the asymptotic series stops only when a term
     exp(-2 k^2 t^2) drops below 1e-10. At t = 0 every term is 1, so the
     loop never ends. At the t of about 1e-16 that cause 1 leaves instead
     of 0 it needs more than 1e16 terms, and the term counter `k` is a
     double that stops growing at 2^53 (9.0e15), so it never ends either
     (the grid's two such points need about 2.2e16 terms); with larger
     noise (about 1e-15) it may end after weeks. Without ties D is never 0, so
     only samples with repeated values reach this, for example y = x.
  3. **Exact p-value as 1 - cdf** (10 points, n1 = 8 to 100, relative
     errors from 5.4e-11 to 2.7e44): `TwoSidedExact` returns
     `1 - valid_paths / total_paths`, which keeps only an absolute accuracy
     of about 1e-14 (the path counts go up to 9e58 and are added in doubles).
     Every p-value below about 1e-3 loses digits; below about 1e-14 the
     result is rounding noise: 2.8e-15 instead of 2.0e-29, and -6.0e-15 (a
     negative probability) instead of 2.2e-59.
  4. **Asymptotic series stopped at an absolute 1e-10** (4 points, n1 = 50
     to 101, up to 5.4e-11 off): the same function as
     [finding 053](053-statrs-ks-onesample.md)'s cause 4. Where
     d sqrt(m n / (m + n)) is small and the true p-value is 1, the result
     is 0.99999999994558 or 1.0000000000293348 (above 1).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::stats_tests::NaNPolicy;
use statrs::stats_tests::ks_test::{KSTwoSampleAlternativeMethod as M, ks_twosample};
use std::sync::mpsc;
use std::time::Duration;

fn ks(x: Vec<f64>, y: Vec<f64>, m: M) {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || tx.send(ks_twosample(x, y, m, NaNPolicy::Error)));
    match rx.recv_timeout(Duration::from_secs(5)) {
        Ok(r) => println!("{:?}", r.unwrap()),
        Err(_) => println!("no result after 5 s"),
    }
}

fn main() {
    let x = vec![0.0, 800.0, 800.0];
    ks(x.clone(), [x.clone(), x.clone()].concat(), M::GreaterAsymptotic);
    ks(x.clone(), x.clone(), M::TwoSidedAsymptotic);
    let low: Vec<f64> = (0..100).map(f64::from).collect();
    let high: Vec<f64> = (1000..1100).map(f64::from).collect();
    ks(low, high, M::TwoSidedExact);
    let a: Vec<f64> = (1..=50).map(|i| f64::from(100 * i)).collect();
    let b: Vec<f64> = (1..=50).map(|j| f64::from(100 * j - 1)).collect();
    ks(a, b, M::TwoSidedAsymptotic);
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below (statistic, p-value). The first,
second and fourth rows are grid points; the third has the same sizes and
statistic as two grid points, so the same p-value.

## Expected

| Samples, method | Cause | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|
| x = (0, 800, 800), y = x twice, greater | 1 | 1.1102230246251565e-16, 0.9999999999999999 | 0, 1 | 0, 1 |
| x = y = (0, 800, 800), two-sided asymptotic | 2 | no result after 5 s | 0, 1 | 0, 1 |
| x = 0..99, y = 1000..1099, two-sided exact | 3 | 1.0000000000000007, -5.995204332975845e-15 | 1, 2.2087606931995028e-59 | 1, 2.2087606931995054e-59 |
| x_i = 100 i, y_j = 100 j - 1 (50 each), two-sided asymptotic | 4 | 0.020000000000000018, 0.999999999963176 | 0.02, 1 | 0.020000000000000018, 1 |

Second opinion: SciPy (`scipy.stats.ks_2samp` with the matching
`alternative` and `method`; for method 3 `scipy.stats.kstwobign.sf` at
SciPy's statistic, because SciPy's own two-sided asymptotic p-value uses a
different approximation) agrees with mpmath within the tolerance at every
point of the dataset except 3 tied points, where SciPy's statistic is
1.1e-16 instead of 0: the report counts 2 of them as "both wrong" (one
cause 1 and one cause 2 point) and the third as "SciPy wrong", where
statrs gives exactly 0. Its worst error elsewhere is 1.4e-14. The exact p-value
also agrees with the Gnedenko and Korolyuk closed form for equal sizes and
with enumeration of every order for small sizes (reference tests).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The routes: the statistic
   exactly in rational arithmetic; the exact p-value by counting the lattice
   paths that stay inside the band in exact integers; Hodges's formula and
   the limit distribution (the ks_1samp series, 1 at t = 0) at the working
   precision. New tests in `reference/tests/test_functions.py` check the
   exact p-value against enumeration of all orders (up to 9 values in all)
   and the Gnedenko and Korolyuk sum for n1 = n2 up to 100 (exact equality
   of rationals), the statistic against walking the pooled sample, D = 0
   for the repeated samples, Hodges's formula against the exact one-sided
   tail (a sanity check of the m, n roles: within 25%, and closer than with
   the sizes swapped), and the rounding bound against the tolerance reason.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random pairs (NumPy `default_rng(20261102)`), n1 and
   n2 log-uniform in 1..100: 70% distinct normal samples (y shifted and
   scaled at random) with any method, 20% integer samples with ties (0 to
   at most 10) and 10% y = x repeated 1 to 5 times, both with methods 0, 1
   and 3. Each point judged at the larger of 5e-12 and its own rounding
   bound (`ks_2samp_rounding_bound`), true values from the reference at 30
   and 50 digits (which agree at every point). The method 3 calls ran in a
   copy of statrs's series with an iteration counter (identical results to
   statrs wherever it returns). 330 fail:
   - 193 through cause 1 (methods 0 and 1, true D = 0, statrs 1.4e-17 to
     7.7e-15); no statistic with a true value above 0 fails at these sizes
     (up to 100). Larger samples (12 normal pairs, `default_rng(20261103)`,
     sizes (1000, 1000) to (100000, 30000), method 3) show that cause 1 also
     moves statistics above 0 once the sums get long: up to 4.0e-11 relative
     at (30000, 100000), beyond 5e-12 and beyond that point's own rounding
     bound (1.5e-11); at (1000, 1000) the error is below 1e-15.
   - 102 through cause 2 (method 3, all with true D = 0: 39 where statrs's
     D is exactly 0, 63 where it is 1.1e-16 to 1.3e-15).
   - 28 through cause 3 (true p from 1.1e-10 to 6.5e-3, absolute errors up
     to 5.1e-14).
   - 7 through cause 4 (results from 0.99999999993 to 1.00000000006, true 1).
5. **Root cause** (from the source, `src/stats_tests/ks_test.rs`, confirmed
   with a copy that has the remedies below; the copy passes every grid
   point and every random point):
   - **Cause 1:** `f1 += 1.0 / n1` and `f2 += 1.0 / n2` per value. Setting
     `f1 = i / n1` and `f2 = j / n2` from the counts removes the cause: equal
     rationals give equal doubles, so D = 0 exactly where it should be.
   - **Cause 2:** `if term.abs() < 1e-10 { break; }` in
     `onesample_kolmogorov_twosided_pvalue` with no other exit: at x = 0 the
     term is always 1. Returning 1 for x = 0 (together with cause 1's
     remedy, which makes D exactly 0) and cause 4's remedy remove it.
   - **Cause 3:** `1.0 - valid_paths / total_paths` in
     `twosample_schroer_and_trenkler_twosided_pvalue`. Counting the paths
     that leave the band instead (a point outside the band starts with all
     C(x + y, x) paths to it, every other point adds its two neighbours) and
     dividing that count by the total gives P(D >= d) without a
     subtraction; SciPy counts the outside paths in the same spirit.
   - **Cause 4:** as finding 053 (stop at a term below 2^-52 times the sum,
     and use the theta form for t < 1).

## Notes

Cause 1 also moves the statistic by a few units in the last place where it
is not 0 (1.0000000000000007 instead of 1 at n1 = n2 = 100), far below the
tolerance; those points pass.

statrs documents `TwoSidedExact` as using Schröer and Trenkler (1995):
"This paper builds on top of Hodges (1957) accounting for ties." The code
does not use ties: the exact routine receives only d and the two sizes, so
it counts the orders of an untied pooled sample (as SciPy does), and unlike
`ks_onesample` it does not refuse tied data. With ties its p-value is
therefore the untied one. Tied samples are not run with method 2 on the
grid, so this is not judged; it is seen in the source only.

Unlike `ks_onesample`, `ks_twosample` never clamps the p-value to [0, 1],
which is why -6.0e-15 (cause 3) and 1.0000000000293348 (cause 4) come out.

Hodges's eq. 5.3 and the limit distribution are approximations by design;
the functions checked are those formulas, not the finite-sample
distribution, so their accuracy as approximations is not judged (it is
within 25% at the sizes in the sanity test). statrs returns
`ExactAndTooLarge` for n1 n2 > 10000 (documented); the grid's (101, 100)
samples skip the exact method.
