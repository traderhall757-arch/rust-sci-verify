# 015: special Beta::inv_inc_beta loops forever, panics, and gives wrong small quantiles

- **Crate / version:** special 0.14.2 (`special::Beta::inv_inc_beta`, the
  inverse of the regularized incomplete beta function: the x with
  I_x(a, b) = p; default features)
- **Severity:** high (an endless loop; silently wrong answers, off by up to
  a factor of 1e188 for small quantiles; a panic in debug builds)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 462 of 1620 points of
  `reference/data/betaincinv.jsonl` (tolerance 2e-13); the exact list is
  `INV_INC_BETA` in `crates/check-special/src/expected.rs`, the report is
  [reports/special-0.14.2.md](../reports/special-0.14.2.md). By cause
  (section "Root cause"):
  1. **p below about 1.5e-162** (161 points: p = 5e-324 and 1e-300 for every
     a, b; the 162nd, a = 1, b = 2, p = 5e-324, is tagged
     `unstable-reference`): in a **release** build the call never returns;
     in a **debug** build (as in CI) it panics with "assertion failed:
     x >= 0.0 && x <= 1.0 && p > 0.0 && q > 0.0" (a `debug_assert!` in
     `inc_beta`). The harness records the panics.
  2. **Endless loop in the far tails** (7 points, both builds): p = 1e-100
     with (a, b) = (100, 2), (1000, 0.1), (1000, 2), (1000, 5), (1e4, 0.1),
     (1e4, 2), (1e4, 5). The harness records "no result within 5 s"
     (`Check::run_with_time_limit`); checked for 60 s.
  3. **Small quantiles wrong** (109 points, off by more than 1e-3): every
     true quantile here is below 4.6e-14. For 9 of them it is below the
     smallest double (0) and the crate returns about 5e-17; for the other
     100 the crate is off by 1.2e-3 up to a factor of 1.7e188 (1.3e-16
     instead of 7.9e-205 at a = 0.5, b = 1e4, p = 1e-100).
  4. **Smaller quantiles inaccurate** (16 points, 1.2e-8 to 1.9e-4): true
     quantiles from 6.1e-15 to 3.3e-4, e.g. 1.9e-4 at a = 0.5, b = 1e4,
     p = 1e-5 (x = 7.85e-15).
  5. **Small errors** (169 points, 2.2e-13 to 5.8e-9): mostly
     [finding 014](014-special-inc-beta.md)'s error in `inc_beta`, which
     the Newton iteration inverts, plus the stopping rule of item 4.

`inc_beta` has its own accuracy problem, [finding 014](014-special-inc-beta.md).
statrs's `inv_beta_reg` is a port of this function and has the same
endless loop, start-value failure and stopping rule
([finding 011](011-statrs-students-t-inverse-cdf.md)); this finding is the
original.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add special@=0.14.2
use special::Beta;
use std::{sync::mpsc, thread, time::Duration};

/// The x with I_x(a, b) = p, with ln_beta precomputed as the documentation asks.
/// Runs on its own thread: some calls never return.
fn inv_inc_beta(a: f64, b: f64, p: f64) -> Result<f64, &'static str> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || tx.send(p.inv_inc_beta(a, b, a.ln_beta(b))));
    match rx.recv_timeout(Duration::from_secs(60)) {
        Ok(x) => Ok(x),
        Err(mpsc::RecvTimeoutError::Timeout) => Err("no result within 60 s"),
        Err(mpsc::RecvTimeoutError::Disconnected) => Err("panicked"),
    }
}

fn main() {
    println!("{:?}", inv_inc_beta(2.0, 1.0, 1e-300)); // should be 1e-150 (= p^(1/a))
    println!("{:?}", inv_inc_beta(1000.0, 2.0, 1e-100)); // should be 0.7900887907081466
    println!("{:?}", inv_inc_beta(0.5, 0.1, 1e-10)); // should be 3.2053074662075163e-19
    println!("{:?}", inv_inc_beta(0.1, 0.1, 1e-100)); // should be 0 (true value about 1e-997)
    println!("{:?}", inv_inc_beta(0.5, 1e4, 1e-5)); // should be 7.854177986380854e-15
    println!("{:?}", inv_inc_beta(2.0, 1e4, 0.01)); // should be 1.485462097727329e-5
    println!("{:?}", inv_inc_beta(1e4, 0.1, 0.5)); // should be 0.999999940658221
}
```

This is the program that was run. In release mode it prints
`Err("no result within 60 s")` twice and then the values of the "Crate"
column below; in debug mode the first line is `Err("panicked")` (after the
assertion message) and the rest is the same.

## Expected

| a | b | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 2 | 1 | 1e-300 | never returns (release); panic (debug) | 1e-150 | 1e-150 | no result |
| 1000 | 2 | 1e-100 | never returns | 0.7900887907081466 | 0.7900887907081466 | no result |
| 0.5 | 0.1 | 1e-10 | 1.3528601681509085e-16 | 3.2053074662075163e-19 | 3.205307466207515e-19 | 421 |
| 0.1 | 0.1 | 1e-100 | 5.5494638638101455e-17 | 0 (about 1e-997) | 2.225073858507201e-308 | 2.5e291 (relative to 2.2e-308) |
| 0.5 | 1e4 | 1e-5 | 7.852715030242036e-15 | 7.854177986380854e-15 | 7.854177986380854e-15 | 1.9e-4 |
| 2 | 1e4 | 0.01 | 1.4854621063997174e-5 | 1.485462097727329e-5 | 1.4854620977273294e-5 | 5.8e-9 |
| 1e4 | 0.1 | 0.5 | 0.9999999457720979 | 0.999999940658221 | 0.999999940658221 | 5.1e-9 |

Second opinion: SciPy (`scipy.special.betaincinv`) agrees with mpmath within
4.6e-16 at every row except (0.1, 0.1, 1e-100). There the lower-tail power
law I_x(a, b) ~ x^a / (a B(a, b)) gives x ~ (p a B(a, b))^(1/a) =
(1e-100 * 0.1 * 19.71)^10, about 1e-997, so the true quantile is 0 as a
double; SciPy's 2.2e-308 is also wrong (its tolerance reason lists such
points). At (2, 1, 1e-300) the closed form x = p^(1/a) (tested in
`reference/tests/test_functions.py`) gives 1e-150.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `special = "=0.14.2"`
   (default features, as in the check crate; the program above, in release
   and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The quantile is Newton's method on log I (or
   log(1 - I)) in s = log(x / (1 - x)) after bisection, with I from the
   continued fraction DLMF 8.17.22; it is tested by I(x) = p with I from
   the independent positive-term series (including the rows above where
   SciPy is wrong) and against closed forms for b = 1, a = 1 and
   a = b = 1/2 in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261014)`):
   a, b = 10^u with u uniform in [-1, 5]; p = 10^v with v uniform in
   [-320, log10(1/2)] for half the points and p uniform in (0, 1) for the
   other half; each rounded to 6 significant digits; mpmath at 50 and 100
   digits (they agreed at all 3,000). Each call ran in its own process
   (release build) with a 2 s limit. 1,520 fail (51%):
   - 809 never return: exactly the 763 points with min(p, 1 - p) below
     1.5e-162, and 46 more in the far lower tail (p from 2.7e-40 down,
     max(a, b) >= 130);
   - 204 off by more than 1e-3, all with true quantiles below 1.6e-7;
   - 50 off by 1e-8 to 7.2e-4 (true quantiles up to 1.6e-3; e.g. 7.2e-4 at
     a = 0.204157, b = 71836.3, p = 0.752562, x = 2.6e-6);
   - 457 off by 2e-13 to 8.3e-9.
5. **Root cause** (confirmed with a copy of `inv_inc_beta`, scratch, not
   committed, that reproduces the crate's results and failure counts on the
   whole grid and can switch the pieces below). `inv_inc_beta` is AS 64 /
   AS 109 with the AS R83 stopping rule: a starting value, then Newton
   steps (I(x) - p) / f(x) whose size is cut by factors of 3 until the
   step's square is below the previous one and stays inside [0, 1].
   - **Start value NaN** (item 1): the start uses
     `(-(a * a).ln()).sqrt()` with a = p (or 1 - p for p > 1/2); a * a
     underflows to 0 for a below about 2^-537.5 = 1.5e-162, so the start is
     NaN. The clamp to [1e-4, 0.9999] uses `<` comparisons, which are false
     for NaN, so `inc_beta(NaN)` is called: its `debug_assert!` panics in a
     debug build, and in a release build its series never meets its
     stopping test (comparisons with NaN are false), so it loops forever.
     Computing the start as `(-2.0 * a.ln()).sqrt()` removes all 162 grid
     failures of this kind.
   - **Endless step loop** (item 2): the Newton step divides by the beta
     density, written as `exp(ln_beta + (1 - p) ln x + (1 - q) ln(1 - x))`.
     Far in the tail of a large parameter, I(x) at the start value
     underflows to 0 and the exponent exceeds 709 (at a = 1e4, b = 5,
     p = 1e-100 the start is x = 0.836 and the step is -inf). The inner loop
     then divides its factor `g` by 3, but an infinite step never becomes
     small enough; once `g` underflows to 0 the step is NaN. There is no
     iteration limit. With the start fixed, 48 grid points end in this
     loop.
   - **Absolute stopping rule** (items 3 and 4): the iteration stops when
     the step's square is below `acu = 10^e` with
     e = trunc(-13 - 5/p^2 - 1/a^0.2), at least 1e-30 (AS R83): an absolute
     test, so steps below 1e-15 (or, for larger e, 1e-7) are not taken.
     The start value is also clamped to at least 1e-4, and for a small
     parameter Newton's method in x shrinks x by only a constant factor per
     step (I ~ x^a is far from linear), so tiny quantiles are not reached.
     At a = 0.1, b = 0.1, p = 1e-100 the copy stops at 5.5e-17. With the
     start fixed and a relative stopping test instead (at most 200 steps),
     passes rise from 1157 to 1275 of 1620, but 154 points are still off
     by more than 1e-3: both the absolute test and the iteration in x
     itself limit small quantiles.
   - **`inc_beta`'s error** (item 5): Newton's method can only be as
     accurate as I(x). For 131 of the 143 small errors that remain in the
     copy with the start fixed and the relative stopping test, the error
     equals finding 014's error of `inc_beta` at the true quantile divided
     by the condition number x f(x) / I (within a factor of 5, mostly to
     2 or 3 digits).
