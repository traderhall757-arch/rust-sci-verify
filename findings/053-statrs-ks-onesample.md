# 053: statrs ks_onesample gives wrong p-values: an exact tail taken as 1 - cdf, n^n and binomial overflow, and an asymptotic series stopped at 1e-10

- **Crate / version:** statrs 0.19.1 (`statrs::stats_tests::ks_test::ks_onesample`,
  all five `KSOneSampleAlternativeMethod`s)
- **Severity:** high (silently wrong answers: a two-sided exact p-value of
  2.7e-15 instead of 6.6e-27; NaN or 1 for samples of 144 to 169 values; a
  one-sided p-value of 1 instead of 0.117 for 1100 values)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 30 of 544 points of `reference/data/ks_1samp.jsonl`
  (tolerance 5e-12, for the statistic and the p-value; samples tested
  against the uniform distribution on (0, 1), whose cdf statrs computes
  exactly). The exact list is `KS_1SAMP` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). The statistic is
  always right. Four causes (section "Root cause"):
  1. **Two-sided exact p-value as 1 - cdf** (12 points, n = 20 to 100,
     relative errors from 7.2e-12 to 4.1e11): `TwoSidedExact` returns
     `1 - P(D < d)`, which keeps only an absolute accuracy of about 3e-15.
     Every p-value below about 1e-3 loses digits; below about 1e-14 the
     result is rounding noise (2.7e-15 where the true value is 6.6e-27).
  2. **n^n overflows from n = 144** (11 points, all at n = 150): the exact
     method scales by `factorial(n) / n.powi(n)`, and 144^144 exceeds the
     largest double. The p-value is NaN for most samples of 144 to 169
     values (statrs refuses 170 and more); where n! times the matrix entry
     stays finite (small statistics, seen at n = 144 and 145) it is exactly
     1 instead, also wrong unless the true p-value rounds to 1 (n = 144,
     D = 0.0285: 1 instead of 0.99962; off the grid).
  3. **Binomial coefficients overflow from n = 1030** (3 points, all at
     n = 1100): `Less`, `Greater` and `TwoSidedApproximate` sum the
     Birnbaum and Tingey series in plain doubles, and C(n, j) is infinite
     for j near n / 2 once n >= 1030. The p-value is then 1 or NaN.
     Below that, the terms' factors can underflow before they are
     multiplied by a large C(n, j): from about n = 700, p-values below
     about 1e-89 lose digits or become 0 (seen off the grid only).
  4. **Asymptotic series stopped at an absolute 1e-10** (4 points, n = 10
     to 150, up to 6.1e-11 off): `TwoSidedAsymptotic` stops when a term
     drops below 1e-10. For D sqrt(n) below about 0.25, where the true
     p-value is 1 to double precision, many terms are near 1 and the
     truncation error stays near 1e-10: 0.99999999993886 instead of 1.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::Uniform;
use statrs::stats_tests::NaNPolicy;
use statrs::stats_tests::ks_test::{KSOneSampleAlternativeMethod as M, ks_onesample};

/// z_i = ((7919 (i + 1)) mod 2048 + 1/2) / 2048, times c.
fn sample(n: usize, c: f64) -> Vec<f64> {
    (0..n).map(|i| c * (((7919 * (i + 1)) % 2048) as f64 + 0.5) / 2048.0).collect()
}

fn main() {
    let u = Uniform::new(0.0, 1.0).unwrap();
    let ks = |x: Vec<f64>, m: M| ks_onesample(x, &u, m, NaNPolicy::Error).unwrap();
    println!("{:?}", ks(sample(100, 0.5), M::TwoSidedExact));
    println!("{:?}", ks(sample(150, 1.0), M::TwoSidedExact));
    println!("{:?}", ks(sample(1100, 1.0), M::Greater));
    println!("{:?}", ks(sample(1100, 1.0), M::TwoSidedApproximate));
    let spread: Vec<f64> = (1..=100).map(|i| (2 * i - 1) as f64 / 200.0).collect();
    println!("{:?}", ks(spread, M::TwoSidedAsymptotic));
}
```

This is the program that was run (release and debug mode give the same
output); it gives the "Crate" column below (statistic, p-value).

## Expected

| Sample, method | Cause | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|
| 100 values pushed to 0 (c = 1/2), two-sided exact | 1 | 0.5316162109375, 2.6645352591003757e-15 | 0.5316162109375, 6.560340895851515e-27 | 0.5316162109375, 6.560340895851515e-27 |
| 150 values (c = 1), two-sided exact | 2 | 0.06199544270833335, NaN | 0.061995442708333334, 0.5895570165667271 | 0.06199544270833335, 0.5895558723971797 |
| 1100 values (c = 1), greater | 3 | 0.031048473011363598, 1.0 | 0.031048473011363636, 0.11747311668035498 | 0.031048473011363598, 0.11747311668035561 |
| 1100 values (c = 1), two-sided approximate | 3 | 0.031048473011363598, 1.0 | 0.031048473011363636, 0.23494623336070997 | 0.031048473011363598, 0.23494623336071122 |
| x_i = (2 i - 1) / 200, two-sided asymptotic | 4 | 0.0050000000000000044, 0.999999999938857 | 0.005000000000000053, 1.0 | 0.0050000000000000044, 1.0 |

Second opinion: SciPy (`scipy.stats.ks_1samp` with the matching
`alternative` and `method`) agrees with mpmath within the tolerance at
every point of the dataset except 4 two-sided exact points at n = 150,
including the second row (1.9e-6 off): above n = 140 SciPy's `kstwo`
switches to an approximation. There, Steck's (1971) determinant in exact
rational arithmetic agrees with mpmath to 38 digits at every n = 150 point
(`test_ks_two_sided_exact_matches_steck_where_scipy_is_wrong`). The
statistics differ in the last digits only (1.4e-14 relative at most,
n = 1100, method 0; both round i / n).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double, for both outputs. The routes: the statistic
   exactly in rational arithmetic; the one-sided p-value (and the two-sided
   one for d >= 1/2, twice it) from the Birnbaum and Tingey sum in exact
   integer arithmetic; the two-sided exact one below 1/2 from the
   Marsaglia, Tsang and Wang matrix (no negative entries) with the digits
   that 1 - P loses added first; the limit distribution from its two
   series. New tests in `reference/tests/test_functions.py` check the
   two-sided and one-sided p-values against Steck's determinant in exact
   rational arithmetic (n up to 150), Marsaglia, Tsang and Wang's published
   K(10, 0.274) = 0.6284796154565043, the closed forms for n = 1,
   d >= 1 - 1/n and d <= 1/n, the matrix route against the d >= 1/2
   identity, the limit series against `mpmath.jtheta`, and the statistic
   against counting the empirical cdf.
3. **Second opinion**: see "Expected".
4. **Off-grid**: two scans, each point judged at the larger of 5e-12 and
   its own rounding bound (`ks_1samp_rounding_bound`), with the true values
   from the reference at 30 and 50 digits (which agree at every point):
   - 3,000 random samples (NumPy `default_rng(20261031)`): 80% with n
     log-uniform in 1..169 and any of the five methods, 20% with n uniform
     in 170..2999 and the four other methods; values u^g with u uniform on
     (0, 1) and g = exp(N(0, 0.6)). 372 fail: 45 through cause 1 (two-sided
     exact, n < 144, true p from 1.4e-33 to 2.3e-4; statrs's absolute error
     is at most 2.7e-15 at every exact point below n = 144), 19 through
     cause 2 (n = 148 to 169, all NaN), 301 through cause 3 (n from 1043,
     NaN at 222, 1 at 79) and 7 through cause 3's underflow (n = 706 to
     1026, true p from 2.0e-89 to 1.0e-203, relative errors 2.1e-9 to 1).
   - 1,000 nearly even samples (`default_rng(20261101)`): x_i =
     (i - 1/2 + u_i) / n, u_i uniform on (-w, w), w uniform on (0, 1/2),
     n uniform in 2..169, any method. 93 fail, all two-sided asymptotic
     (cause 4), at D sqrt(n) from 0.041 to 0.24, up to 6.8e-11 off.
5. **Root cause** (from the source, `src/stats_tests/ks_test.rs`, confirmed
   with a copy that has the remedies below):
   - **Cause 1:** `1.0 - onesample_marsaglia_et_al_twosided_pvalue(...)`.
     The matrix result P(D < d) itself is accurate (the absolute error
     of 1 - P stays at most 4.0e-15 on both the grid and the random scan),
     so the loss is the subtraction. For d >= 1/2 the one-sided events are
     disjoint and P(D >= d) = 2 P(D+ >= d) exactly, which removes the
     cause there (SciPy does the same); below 1/2 the upper tail needs a
     method of its own (Simard and L'Ecuyer 2011 discuss the choices). The
     copy does not have one, so 11 grid points still fail in it through
     cause 1 (the 7 cause 1 points with d < 1/2, and 4 cause 2 points at
     n = 150 with d < 1/2 and p below 1e-8) and 40 random ones, all with
     absolute errors at most 4.0e-15.
   - **Cause 2:** `v[k] * factorial::factorial(n as u64) / n.powi(n as i32)`:
     n.powi(n) is infinite from n = 144 (143^143 = 1.6e308), and
     v[k] n! is infinite too, so the result is inf / inf = NaN; where
     v[k] n! stays finite it is finite / inf = 0, and the p-value 1. Forming
     n! / n^n as the product of i / n (at least 1e-74 for n <= 169) removes
     the cause: no NaN is left on the grid or in the scan, and the
     n = 144 and 145 samples above give the true p-value.
   - **Cause 3:** `factorial::binomial(n as u64, j)` is infinite once
     C(n, j) > 1.8e308 (from C(1030, 515)); times a positive power it
     gives inf, times an underflowed 0 it gives NaN, and the final clamp
     maps inf to 1 and leaves NaN. Below n = 1030, the factors
     `(j / n + d).powi(j - 1)` and `(1 - d - j / n).powi(n - j)` can
     underflow although their product with C(n, j) is a normal number.
     Carrying an exponent alongside each factor (or summing the terms'
     logarithms with a scaled sum) removes both.
   - **Cause 4:** `if term.abs() < 1e-10 { break; }` in
     `onesample_kolmogorov_twosided_pvalue`. Stopping at a term below
     2^-52 times the sum, and using the theta form 1 - sqrt(2 pi) / t
     sum exp(-(2k - 1)^2 pi^2 / (8 t^2)) for t < 1 (where the alternating
     series has many terms), removes the cause.

   The copy with the remedies for causes 2 to 4 passes every grid point
   and every random point apart from the cause 1 points below d = 1/2 and
   2 random points with true p below 1e-303, where its own terms are
   subnormal (statrs gives NaN at both: n = 1713 and 1819, cause 3). The
   nearly even scan passes completely.

## Notes

Cause 1 also moves points that pass: at n = 50, x = (3/4) z, two-sided
exact, statrs is 4.2e-12 off at p = 1.86e-4, 37 times that point's own
rounding bound (1.1e-13), but below the dataset's single tolerance of
5e-12, so it is not listed.

statrs documents `TwoSidedExact` as producing "the exact p-value" for
n < 140 and `TwoSidedAsymptotic` as having "roughly 5 digits of precision"
for n > 140. The limit distribution is the function checked here, not the
finite-n distribution, so its accuracy as an approximation is not judged;
only the computation of 1 - K(d sqrt(n)) is, and that is right to about
1e-10 (cause 4), not to double precision.

Not on the grid (no listed points): samples with ties make `TwoSidedExact`
return an error (documented), and samples of 170 or more values make it
return `ExactAndTooLarge` (documented). Only the uniform distribution is
used; with another distribution the statistic also carries that
distribution's cdf error (for example finding 009 for `Normal`).
