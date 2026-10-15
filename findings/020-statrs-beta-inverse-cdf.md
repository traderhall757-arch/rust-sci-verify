# 020: statrs Beta inverse_cdf panics, loops forever, and gives wrong small quantiles

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Beta`, method
  `inverse_cdf` of `ContinuousCDF`, which calls
  `statrs::function::beta::inv_beta_reg`)
- **Severity:** high (an endless loop; silently wrong answers, off by up to
  a factor of 1.7e188 for small quantiles; a panic for every p below about
  1.5e-162)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 479 of 1620 points of
  `reference/data/betaincinv.jsonl` (tolerance 2e-13); the exact list is
  `BETA_PPF` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). By cause
  (section "Root cause"):
  1. **p below about 1.5e-162** (161 points: p = 5e-324 and 1e-300 for every
     a, b; the 162nd, a = 1, b = 2, p = 5e-324, is tagged
     `unstable-reference`): the call **panics** with "called
     `Result::unwrap()` on an `Err` value: XOutOfRange", in release and
     debug builds alike.
  2. **Endless loop in the far tails** (7 points): p = 1e-100 with
     (a, b) = (100, 2), (1000, 0.1), (1000, 2), (1000, 5), (1e4, 0.1),
     (1e4, 2), (1e4, 5). The harness records "no result within 5 s"
     (`Check::run_with_time_limit`); checked for 60 s.
  3. **Small quantiles wrong** (109 points, off by more than 1e-3): every
     true quantile here is below 4.6e-14. For 9 of them it is below the
     smallest double (0) and the crate returns about 5e-17; for the other
     100 the crate is off by 1.2e-3 up to a factor of 1.7e188 (1.3e-16
     instead of 7.9e-205 at a = 0.5, b = 1e4, p = 1e-100).
  4. **Smaller quantiles inaccurate** (16 points, 1.1e-8 to 1.9e-4): true
     quantiles from 6.1e-15 to 3.3e-4, e.g. 1.9e-4 at a = 0.5, b = 1e4,
     p = 1e-5 (x = 7.85e-15).
  5. **Small errors** (186 points, 2.0e-13 to 5.8e-9): mostly
     [finding 019](019-statrs-beta-cdf.md)'s error in `beta_reg`, which
     the Newton iteration inverts, plus the stopping rule of item 4.

`inv_beta_reg` is a port of special's `inv_inc_beta` (its source says so)
and has the same counts for items 1 to 4 as
[finding 015](015-special-inv-inc-beta.md), which has the detailed
analysis; all 293 such points are also in that finding's list. Only the
symptom of item 1 differs (special loops forever in a release build). `StudentsT::inverse_cdf` uses the same function
([finding 011](011-statrs-students-t-inverse-cdf.md)).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Beta, ContinuousCDF};
use std::{sync::mpsc, thread, time::Duration};

/// The x with I_x(a, b) = p. Runs on its own thread: some calls never return.
fn inverse_cdf(a: f64, b: f64, p: f64) -> Result<f64, &'static str> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || tx.send(Beta::new(a, b).unwrap().inverse_cdf(p)));
    match rx.recv_timeout(Duration::from_secs(60)) {
        Ok(x) => Ok(x),
        Err(mpsc::RecvTimeoutError::Timeout) => Err("no result within 60 s"),
        Err(mpsc::RecvTimeoutError::Disconnected) => Err("panicked"),
    }
}

fn main() {
    println!("{:?}", inverse_cdf(2.0, 1.0, 1e-300)); // should be 1e-150 (= p^(1/a))
    println!("{:?}", inverse_cdf(1000.0, 2.0, 1e-100)); // should be 0.7900887907081466
    println!("{:?}", inverse_cdf(0.5, 0.1, 1e-10)); // should be 3.2053074662075163e-19
    println!("{:?}", inverse_cdf(0.1, 0.1, 1e-100)); // should be 0 (true value about 1e-997)
    println!("{:?}", inverse_cdf(0.5, 1e4, 1e-5)); // should be 7.854177986380854e-15
    println!("{:?}", inverse_cdf(2.0, 1e4, 0.01)); // should be 1.485462097727329e-5
    println!("{:?}", inverse_cdf(1e4, 0.1, 0.5)); // should be 0.999999940658221
}
```

This is the program that was run. In release and debug mode alike it
prints the panic message and `Err("panicked")`, then
`Err("no result within 60 s")`, then the values of the "Crate" column
below.

## Expected

| a | b | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 2 | 1 | 1e-300 | panic | 1e-150 | 1e-150 | no result |
| 1000 | 2 | 1e-100 | never returns | 0.7900887907081466 | 0.7900887907081466 | no result |
| 0.5 | 0.1 | 1e-10 | 1.3528601681509346e-16 | 3.2053074662075163e-19 | 3.205307466207515e-19 | 421 |
| 0.1 | 0.1 | 1e-100 | 5.549463863810121e-17 | 0 (about 1e-997) | 2.225073858507201e-308 | 2.5e291 (relative to 2.2e-308) |
| 0.5 | 1e4 | 1e-5 | 7.852715029998704e-15 | 7.854177986380854e-15 | 7.854177986380854e-15 | 1.9e-4 |
| 2 | 1e4 | 0.01 | 1.485462106399717e-5 | 1.485462097727329e-5 | 1.4854620977273294e-5 | 5.8e-9 |
| 1e4 | 0.1 | 0.5 | 0.9999999457720979 | 0.999999940658221 | 0.999999940658221 | 5.1e-9 |

Second opinion: as in finding 015 (the same rows): SciPy
(`scipy.special.betaincinv`) agrees with mpmath within 4.6e-16 at every row
except (0.1, 0.1, 1e-100), where the lower-tail power law gives a quantile
of about 1e-997, 0 as a double, and SciPy's 2.2e-308 is also wrong; at
(2, 1, 1e-300) the closed form x = p^(1/a) (tested in
`reference/tests/test_functions.py`) gives 1e-150.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the `betaincinv` dataset of finding 015: mpmath at
   50 and 100 digits round to the same double; the quantile (Newton's
   method on log I in s = log(x / (1 - x)) with the continued fraction) is
   tested by I(x) = p with the independent positive-term series and against
   closed forms for b = 1, a = 1 and a = b = 1/2 in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261015)`), drawn
   as in finding 015: a, b = 10^u with u uniform in [-1, 5]; p = 10^v with v
   uniform in [-320, log10(1/2)] for half the points (5e-324 where that
   rounds to 0) and p uniform in (0, 1) for the other half; each rounded to
   6 significant digits; mpmath at 50 and 100 digits (they agreed at all
   3,000). Each call ran in its own process (release build) with a 2 s
   limit. 1,552 fail (52%):
   - 748 panic: exactly the 748 points with min(p, 1 - p) below
     2^-537.5 = 1.5e-162;
   - 56 never return, all with p <= 2.0e-30 and max(a, b) >= 88.6;
   - 220 off by more than 1e-3, all with true quantiles below 1.6e-7;
   - 48 off by 1e-8 to 7.0e-4 (e.g. 7.0e-4 at a = 11.0136, b = 27.0904,
     p = 1.6626e-136, x = 7.2e-14);
   - 480 off by 2e-13 to 8.7e-9.
5. **Root cause** (confirmed with a copy of `inv_beta_reg`, scratch, not
   committed, that calls statrs's own `beta_reg`, reproduces
   `Beta::inverse_cdf` bit for bit at all 1,619 untagged grid points and
   can switch the pieces below). `inv_beta_reg` is AS 64 / AS 109 with the
   AS R83 stopping rule, as in finding 015, whose analysis applies line for
   line:
   - **Start value NaN** (item 1): the start uses
     `(-(x * x).ln()).sqrt()`, and x * x underflows to 0 for x below
     2^-537.5, so the start is NaN; `p.clamp(0.0001, 0.9999)` keeps the NaN,
     and `beta_reg(a, b, NaN)` panics because it unwraps
     `checked_beta_reg`'s `XOutOfRange` error (special's `inc_beta` loops
     instead). With the start computed as `(-2.0 * x.ln()).sqrt()` the copy
     has no panics left, on the grid or off it.
   - **Endless step loop** (item 2): the Newton step multiplies by the
     reciprocal of the beta density,
     `exp(ln_beta + (1 - a) ln p + (1 - b) ln(1 - p))` at the current
     iterate p, which overflows far in the tail of a large parameter; the infinite step is
     divided by 3 until its factor underflows, the step becomes NaN, and
     there is no iteration limit. With the start fixed, 48 grid points (and
     192 random points) end in this loop, as in finding 015.
   - **Absolute stopping rule** (items 3 and 4, and the near-1 cases of
     item 5): the iteration stops when the squared step or the squared
     residual is below `acu`, at least 1e-30, an absolute test; the start is
     clamped to at least 1e-4, and Newton's method in x reaches tiny
     quantiles only slowly. With the start fixed and a relative stopping
     test (at most 200 steps) the copy still has 141 grid points off by more
     than 1e-3, so both the stopping rule and the iteration in x limit
     small quantiles.
   - **`beta_reg`'s error** (item 5): for 113 of the 186 small errors the
     error is within a factor of 5 of finding 019's error of `cdf` (or
     `sf` for p > 1/2) at the true quantile divided by the condition number
     x f(x) / p. The others are quantiles near 0 or 1 with b or a small
     (e.g. 5.1e-9 at a = 1e4, b = 0.1, p = 0.5, x = 0.99999994), where the
     absolute stopping rule ends the iteration early.
