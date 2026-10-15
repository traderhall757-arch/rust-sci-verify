# 022: statrs FisherSnedecor inverse_cdf panics, loops forever, and loses the upper quantiles through z / (1 - z)

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::FisherSnedecor`,
  method `inverse_cdf` of `ContinuousCDF`, which calls
  `statrs::function::beta::inv_beta_reg`)
- **Severity:** high (an endless loop; a panic for every p below about
  1.5e-162; silently wrong answers: +inf for a finite quantile, 4.5e15 for
  8.5e38, 0.32% off at the median for d1 = 1e5, d2 = 0.5)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 724 of 1620 points of `reference/data/f_ppf.jsonl`
  (tolerance 2e-13); the exact list is `F_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). By cause
  (section "Root cause"):
  1. **p below about 1.5e-162** (162 points: p = 5e-324 and 1e-300 for every
     d1, d2): the call **panics** with "called `Result::unwrap()` on an
     `Err` value: XOutOfRange", in release and debug builds alike.
  2. **Endless loop** (7 points): p = 1e-100 with (d1, d2) = (1000, 0.5),
     (1000, 1), (1000, 5), (1000, 10), (1e5, 0.5), (1e5, 5), (1e5, 10). The
     harness records "no result within 5 s" (`Check::run_with_time_limit`);
     checked for 60 s.
  3. **Wrong beta quantile** (345 points): `inv_beta_reg`'s own errors
     ([finding 020](020-statrs-beta-inverse-cdf.md), items 3 to 5), carried
     over: 93 points off by more than 1e-3, all with a true beta quantile z
     below 1.6e-15 (e.g. 1.56e-11 instead of 0 at d1 = 0.5, d2 = 1e5,
     p = 1e-100, where the true quantile is 2.7e-400),
     39 from 1e-9 to 4.5e-4 (2.7190e-12 instead of 2.7202e-12 at
     d1 = 0.5, d2 = 100, p = 0.001), and 213 below 1e-9, at all p.
  4. **Upper quantiles lost through z / (1 - z)** (210 points, true
     quantiles from 0.051 to 8.5e58): even the correctly rounded z gives an
     error above the tolerance here. 4 return +inf (z rounded to 1, e.g.
     d1 = 100, d2 = 2, p = 1 - 1e-15, true 1.0e15), 46 are capped near
     2^52 d2 / d1 (4.5e15 instead of 8.5e38 at d1 = d2 = 0.5,
     p = 1 - 1e-10), 67 are off by 1e-9 to 10% and 93 by up to 1e-9. The
     error of `inv_beta_reg`'s z is amplified the same way: 0.32% at the
     median for d1 = 1e5, d2 = 0.5 (5.742344 instead of 5.724210).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, FisherSnedecor};
use std::{sync::mpsc, thread, time::Duration};

/// The x with F(x; d1, d2) = p. Runs on its own thread: some calls never return.
fn inverse_cdf(d1: f64, d2: f64, p: f64) -> Result<f64, &'static str> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || tx.send(FisherSnedecor::new(d1, d2).unwrap().inverse_cdf(p)));
    match rx.recv_timeout(Duration::from_secs(60)) {
        Ok(x) => Ok(x),
        Err(mpsc::RecvTimeoutError::Timeout) => Err("no result within 60 s"),
        Err(mpsc::RecvTimeoutError::Disconnected) => Err("panicked"),
    }
}

fn main() {
    println!("{:?}", inverse_cdf(5.0, 2.0, 1e-300)); // should be 4e-121
    println!("{:?}", inverse_cdf(1000.0, 0.5, 1e-100)); // should be 0.0008794687369934656
    println!("{:?}", inverse_cdf(100.0, 2.0, 0.999999999999999)); // should be 1000799917193443.0
    println!("{:?}", inverse_cdf(0.5, 0.5, 0.9999999999)); // should be 8.462349668255225e38
    println!("{:?}", inverse_cdf(1e5, 0.5, 0.5)); // should be 5.724209992059571
    println!("{:?}", inverse_cdf(1e5, 1e5, 0.5)); // should be 1.0
    println!("{:?}", inverse_cdf(0.5, 1e5, 1e-100)); // should be 0 (true value 2.7e-400)
    println!("{:?}", inverse_cdf(0.5, 100.0, 0.001)); // should be 2.7202381208667023e-12
    println!("{:?}", inverse_cdf(2.0, 1e5, 0.999)); // should be 6.908232471787367
}
```

This is the program that was run. In release and debug mode alike it
prints the panic message and `Err("panicked")`, then
`Err("no result within 60 s")`, then the values of the "Crate" column
below.

## Expected

| d1 | d2 | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 5 | 2 | 1e-300 | panic | 4e-121 | 3.999999999999939e-121 | - |
| 1000 | 0.5 | 1e-100 | no result | 0.0008794687369934656 | 0.0008794687369934658 | - |
| 100 | 2 | 1 - 1e-15 | inf | 1000799917193443.0 | 1000799917193443.0 | inf |
| 0.5 | 0.5 | 1 - 1e-10 | 4503599627370496.0 | 8.462349668255225e38 | 8.462349668255236e38 | 1.0 |
| 1e5 | 0.5 | 0.5 | 5.74234400000132 | 5.724209992059571 | 5.7242099920595715 | 3.2e-3 |
| 1e5 | 1e5 | 0.5 | 0.9999999999893419 | 1.0 | 1.0 | 1.1e-11 |
| 0.5 | 1e5 | 1e-100 | 1.56108214681264e-11 | 0.0 | 4.4501477170144015e-303 | 7.0e296 (relative to 2.2e-308) |
| 0.5 | 100 | 0.001 | 2.7190017822720265e-12 | 2.7202381208667023e-12 | 2.720238120866699e-12 | 4.5e-4 |
| 2 | 1e5 | 0.999 | 6.908232465767259 | 6.908232471787367 | 6.9082324717873655 | 8.7e-10 |

Second opinion: SciPy (`scipy.stats.f.ppf`) agrees with mpmath within
1.6e-14 at every row except the seventh, where it returns 4.5e-303 for a
quantile below the smallest double (one of SciPy's 101 failures on this
grid, all with p <= 1e-100). There, and at the first row, closed forms
confirm mpmath: for d2 = 2 the quantile is (2/d1) z / (1 - z) with
z = p^(2/d1) (4e-121 at row 1), and the seventh row's true value is
2.7e-400 (mpmath at 50 digits; the leading term of the lower tail,
(d2/d1) (p a B(a, b))^(1/a) with a = 1/4, b = 5e4, gives the same).
`reference/tests/test_functions.py` also checks F(ppf(p)) = p through the
independently tested cdf and sf at rows 2, 5 and 7 and at three of SciPy's
failures, and the closed forms for
d1 = d2 = 1, d1 = 2 and d2 = 2 at p from 5e-324 to 1 - 1e-15.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double (they agree at every point). The route
   (`_mp_f_ppf`: x = (d2/d1) e^s with s = log(z / (1 - z)) from Newton's
   method on log I or log(1 - I), `_mp_betaincinv_logit`, the same iteration
   as the `betaincinv` data) is tested as described under "Expected".
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261015)`):
   d1, d2 = 10^u with u uniform in [-0.5, 5]; p uniform in (0, 1) for half
   the points, p = 10^v with v uniform in [-300, 0] for a quarter, and
   p = 1 - 10^v with v uniform in [-16, 0] for the rest; d1, d2 and 10^v
   rounded to 6 significant digits; mpmath at 50 and 100 digits (they agreed
   at all 3,000); statrs in a scratch program with a 5 s limit per call.
   1,615 fail (54%): 329 panics (p from 2.2e-300 to 5.3e-163), 33 calls that
   never return (p from 2.5e-162 to 1.7e-65), 921 errors in z (cause 3; 101
   off by more than 1e-3, all with a true z below 1.2e-6), and 332 of
   cause 4: 330 with z >= 0.9994 (true quantiles from 0.31 up; 26 return
   +inf, 91 more are off by more than 1e-3, the worst 9.1e10 instead of
   4.0e30 at d1 = 46477.5, d2 = 0.941679, p = 1 - 3.1e-15) and 2 with
   subnormal quantiles, where d2 / d1 magnifies the rounding of a subnormal
   z.
5. **Root cause** (from `src/distribution/fisher_snedecor.rs` and
   `src/function/beta.rs`): `inverse_cdf(p)` computes
   `z = inv_beta_reg(d1/2, d2/2, p)` and returns
   `d2 / (d1 * (1/z - 1))`.
   - **Causes 1 to 3** are `inv_beta_reg`'s, analysed in
     [finding 015](015-special-inv-inc-beta.md) (statrs's function is a
     port of special's `inv_inc_beta`) and
     [finding 020](020-statrs-beta-inverse-cdf.md): a NaN start for
     p < 1.5e-162 that `beta_reg` rejects with the unwrapped
     `XOutOfRange` error, a Newton iteration without an iteration limit
     whose step becomes infinite in the far tails, and an absolute stopping
     rule with a start clamped to 1e-4. Classification: at the 345 points of
     cause 3, applying `d2 / (d1 * (1/z - 1))` to the correctly rounded z
     (computed with mpmath) is within tolerance, so the error is in z.
   - **Cause 4, the last step:** the quantile is (d2/d1) z / (1 - z), and
     1/z - 1 = (1 - z) / z is formed from the rounded z. Its relative error
     is about 2^-53 / (1 - z), so it grows with the quantile
     (1 - z = d2 / (d1 x + d2)); it is 2^-52 or 2^-53 once 1 - z is below
     the spacing of doubles near 1 (the cap near 2^52 d2 / d1), and 0 when z
     rounds to 1 (+inf). At the 210 points of cause 4 even the correctly
     rounded z gives an error above the tolerance. Any absolute error in
     `inv_beta_reg`'s z near 1 is amplified the same way, by z / (1 - z)
     (about 1.4e5 at d1 = 1e5, d2 = 0.5, p = 1/2).
     The beta quantile of the other tail, 1 - z = I^-1_(1-p)(d2/2, d1/2),
     computed directly, or Newton's method on the F cdf in log x, avoids
     this.

## Notes

`FisherSnedecor::cdf` and `sf`: [finding 021](021-statrs-fisher-snedecor-cdf.md).
