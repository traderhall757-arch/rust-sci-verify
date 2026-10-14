# 011: statrs StudentsT::inverse_cdf loops forever, panics, and gives wrong quantiles

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::StudentsT`,
  method `inverse_cdf` of `ContinuousCDF`, which calls
  `statrs::function::beta::inv_beta_reg`)
- **Severity:** high (silently wrong answers, off by up to a factor of 1e187
  in the tails and up to 67% in the bulk for large df; an endless loop; a
  panic)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 219 of 513 points of `reference/data/t_ppf.jsonl`
  (tolerance 2e-13); the exact list is `T_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Only location 0
  and scale 1 are checked. By cause (section "Root cause"):
  1. **Endless loop** (23 points): every p of the grid for df = 1e8, and
     (df = 1000, p = 1e-100). The call never returns (checked for 60 s). The
     harness records these as "no result within 5 s" (they are run with
     `Check::run_with_time_limit`). Off the grid it starts at about df = 500
     in the far tails and covers about 73% of all p for df >= 1e7.
  2. **Panic** (38 points): p = 5e-324 and 1e-300 for every df panic with
     "called `Result::unwrap()` on an `Err` value: XOutOfRange". Off the
     grid: every p (or 1 - p) below about 7.4e-163 panics.
  3. **Tails cut off near 1e8** (49 points, plus 2 overflow points): for
     df <= 10, wherever the true |quantile| is above 1e6, the result is
     at most about 1.4e9 in size: e.g. -9.2e7 instead of -3.2e19 (df = 1,
     p = 1e-20) or -8.9e7 instead of -1.6e196 (df = 0.1, p = 1e-20). 42 of
     the 49 are off by more than 50%; the others by 5.9e-11 to 6.4e-2. At the 2 overflow points (df = 0.1,
     p = 1e-100 and 1e-50) the true quantile is below -1.8e308 and statrs
     returns the finite -88896176.33701544.
  4. **Near p = 1/2** (36 points): p = 1/2 +- 1e-10 gives
     +-sqrt(df * 2^-53) (1.05e-8 for df = 1) instead of about 2.5e-10 to
     3.1e-10: up to 4.2e4 times too large (df = 1e6).
  5. **Large df** (51 points, df >= 1e4): errors up to 2.1e-5 (df = 1e5,
     p = 0.25) and 3.1e-3 (df = 1e6, p = 0.25). Off the grid up to 0.67
     (df = 2.75e6, p = 0.067: -0.49995 instead of -1.4974).
  6. **Smaller errors** (20 points): 2.2e-13 to 7.4e-12 at df = 0.5 to 1000
     in the bulk and the tails (e.g. 4.8e-13 at df = 100, p = 0.95).

`StudentsT::cdf` and `sf` have their own problems,
[finding 010](010-statrs-students-t.md); some of the small errors of item 6
may come from `beta_reg` inside the Newton iteration (not separated).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::{sync::mpsc, thread, time::Duration};

fn t(df: f64) -> StudentsT {
    StudentsT::new(0.0, 1.0, df).unwrap()
}

fn main() {
    println!("{:?}", t(1.0).inverse_cdf(1e-20)); // should be -3.1830988618379067e19
    println!("{:?}", t(0.1).inverse_cdf(0.01)); // should be -1.6044257056665652e16
    println!("{:?}", t(0.1).inverse_cdf(1e-50)); // should be -inf (below -1.8e308)
    println!("{:?}", t(1.0).inverse_cdf(0.4999999999)); // should be -3.1415929135263347e-10
    println!("{:?}", t(1e6).inverse_cdf(0.25)); // should be -0.6744899955310873
    println!("{:?}", t(1e5).inverse_cdf(0.05)); // should be -1.6448688647849696
    println!("{:?}", t(100.0).inverse_cdf(0.05)); // should be -1.6602343260853396
    let panicked = std::panic::catch_unwind(|| t(10.0).inverse_cdf(1e-300)).is_err();
    println!("panics: {panicked}"); // should be -2.564525718948198e30
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || tx.send(t(1e8).inverse_cdf(0.25)));
    println!("{:?}", rx.recv_timeout(Duration::from_secs(60))); // should be -0.674489752649431
}
```

This is the program that was run (release mode). It prints the "Crate"
column below, then `panics: true` (with the panic message) and
`Err(Timeout)`.

## Expected

| df | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| 1 | 1e-20 | -92043678.8850747 | -3.1830988618379067e19 | -3.183098861837907e19 | 1.0 |
| 0.1 | 0.01 | -99379281.96647565 | -1.6044257056665652e16 | -1.6044257056665624e16 | 1.0 |
| 0.1 | 1e-50 | -88896176.33701544 | -inf (overflow) | -2.1199605744342963e153 | tagged |
| 1 | 0.4999999999 | -1.0536712127723509e-8 | -3.1415929135263347e-10 | -3.1415942642131064e-10 | 32.5 |
| 1e6 | 0.25 | -0.672383442158132 | -0.6744899955310873 | -0.6744899955310875 | 3.1e-3 |
| 1e5 | 0.05 | -1.644866223383758 | -1.6448688647849696 | -1.6448688647849699 | 1.6e-6 |
| 100 | 0.05 | -1.660234326086132 | -1.6602343260853396 | -1.6602343260853396 | 4.8e-13 |
| 10 | 1e-300 | panic | -2.564525718948198e30 | inf | panic |
| 1e8 | 0.25 | never returns | -0.674489752649431 | -0.6744897526494311 | no result |

Second opinion: SciPy (`scipy.stats.t.ppf`) agrees with mpmath within
1.8e-15 at the rows with df = 0.1, p = 0.01, df = 1, p = 1e-20, and
df >= 100. Where SciPy is also wrong: for df = 1, p = 1/2 - 1e-10 the closed
form tan(pi (p - 1/2)) gives the mpmath value (tested in
`reference/tests/test_functions.py`); for df = 10, p = 1e-300 the mpmath
quantile reproduces p through the cdf to 8e-51; for df = 0.1, p = 1e-50 the
cdf at -1.7976931348623157e308 is 6.2e-32 (by both the continued fraction
and the series DLMF 8.17.8), far above p, so the quantile is below the
largest double.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The quantile is Newton's method on
   log F(-e^u) after bisection, with F as in finding 010; it is tested by
   F(ppf(p)) = p (including df = 1e8 and p = 5e-324) and against the closed
   forms for df = 1 and 2 in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261014)`): df = 10^u with u uniform in [-1.3, 9]; p = 10^v
   with v uniform in [-320, log10(1/2)] for half the points and p uniform in
   (0, 1) for the other half; each rounded to 6 significant digits; mpmath at
   50 and 100 digits (they agreed at all 3,000). Each call ran in its own
   process with a 2 s limit (finding 010's 20,000 points are a separate run
   of the same seed). 2,223 fail (74%): 760 panics (exactly the
   points with min(p, 1 - p) < 7.4e-163), 525 endless loops (from
   df = 498.6 at p = 1.8e-80; 434 of the 598 points with df >= 1e7), 533
   other failures at df >= 1e4 (worst 0.67), 138 tails cut off near 1e8
   (|true quantile| > 1e6), 47 finite values where the true quantile
   overflows, and 220 others (worst 5.5e-6, df = 0.055).
5. **Root cause** (confirmed with an instrumented copy of `inv_beta_reg`,
   scratch, not committed): `inverse_cdf` computes
   y = `inv_beta_reg(df / 2, 0.5, 2 min(p, 1 - p))` and returns
   +-sqrt(df (1 - y) / y). `inv_beta_reg` is AS 64 / AS 109 with the AS R83
   stopping rule (a port from the `special` crate):
   - **Endless loop:** each Newton step is (I_y - target) divided by the beta
     density, written as `exp(ln_beta + (1 - a) ln y + (1 - b) ln(1 - y))`.
     For large a = df/2 the exponent exceeds 709 (4987 at df = 1e8,
     p = 0.25, y = 0.9999; 734.5 at df = 1000, p = 1e-100), so the step is
     -inf. The inner loop then divides its step factor `g` by 3 until the
     step is small enough, which never happens for an infinite step
     (g * q is inf, then NaN once g underflows to 0); there is no iteration
     limit.
   - **Panic:** the starting value uses `(-(x * x).ln()).sqrt()` with
     x = 2p; x * x underflows to 0 for x below about 1.5e-162, so the start
     is NaN, and the first `beta_reg(a, b, NaN)` returns `Err(XOutOfRange)`,
     which `beta_reg` unwraps.
   - **Tails cut off:** the iteration stops on an absolute step size
     (`q * q <= acu` with `acu` at least 1e-30, so steps below about 1e-15
     in y), and the starting value is clamped to [1e-4, 0.9999]. For
     df = 1, p = 1e-20 it stops after 25 iterations at y = 1.18e-16 where
     the answer needs y = 1e-39, giving |t| = sqrt(1 / 1.18e-16) = 9.2e7.
   - **Near p = 1/2:** y is near 1, so 1 - y is 0 or a multiple of 2^-53
     (the double spacing just below 1), and |t| >= sqrt(df * 2^-53): the
     printed -1.0536712127723509e-8 is exactly -sqrt(2^-53) (y = 1 - 2^-53).
   - **Large df:** the AS R83 accuracy is `acu = 10^e` with
     e = trunc(-13 - 5/a^2 - 1/x^0.2) (x = 2 min(p, 1 - p), at least 1e-30),
     here 1e-14, so the iteration stops once the step in y is below
     about 1e-7. With large df
     the answer is in 1 - y, which is about t^2 / df (4.5e-7 at df = 1e6,
     p = 0.25), so an absolute 1e-7 in y is a large relative error in t:
     the copy stops after 6 iterations with 1 - y = 4.521e-7 where
     4.549e-7 is right.

## Notes

Boost.Math computes Student's t quantiles with Hill's algorithm (ACM 396)
for most df and refines with Halley steps on the cdf, and works with both y
and 1 - y, so none of the cancellations above occur. SciPy is also wrong at
36 grid points (report: "both wrong"); its errors differ from statrs's.
