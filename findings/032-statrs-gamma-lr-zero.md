# 032: statrs gamma_lr and gamma_ur panic at x = 0

- **Crate / version:** statrs 0.19.1 (`statrs::function::gamma::gamma_lr`
  and `gamma_ur`, the regularized lower and upper incomplete gamma
  functions; `checked_gamma_lr` / `checked_gamma_ur` return an error)
- **Severity:** medium (a panic instead of the correct value
  P(a, 0) = 0 or Q(a, 0) = 1; no silently wrong value)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 12 of 192 points of `reference/data/gammainc.jsonl`
  (every shape of the grid at x = 0; tolerance 1e-12). The exact list is
  `GAMMA_LR_ZERO` in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Every a > 0
  panics at x = 0, in both functions; x = +inf panics too (finding 024,
  through `InverseGamma`).

`gamma_lr`'s other failures (1 subnormal result on `gammainc` and 14 points
with a >= 1e8 on `gammainc_large`) are finding 012's; the distributions
(`Gamma`, `ChiSquared`, `Erlang`) test x = 0 themselves before calling these
functions, so they are not affected. `Chi` tests x = 0 but then passes
`x * x / 2`, which is 0 for x below 2.72e-162 and +inf from 1.34e154, so it
panics there ([finding 035](035-statrs-chi-cdf.md), added 2026-10-16).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::function::gamma::{checked_gamma_lr, gamma_lr, gamma_ur};

fn main() {
    println!("{:?}", gamma_lr(2.0, 1e-300)); // tiny but fine: 5e-601 rounds to 0.0
    println!("{:?}", checked_gamma_lr(2.0, 0.0)); // Err(XInvalid): P(a, 0) = 0 is rejected
    println!("{:?}", std::panic::catch_unwind(|| gamma_ur(2.0, 0.0)).is_err()); // true: Q(a, 0) = 1 panics
    println!("{:?}", gamma_lr(2.0, 0.0)); // should be 0.0; panics
}
```

This is the program that was run (release and debug mode): it prints
`0.0`, `Err(XInvalid)` and `true` (after the panic message of `gamma_ur`)
and then panics with ``called `Result::unwrap()` on an `Err` value:
XInvalid`` at `src/function/gamma.rs:267`.

## Expected

| a | x | Function | Crate | mpmath (true value) | SciPy |
|---|---|---|---|---|---|
| 2 | 0 | `gamma_lr` | panic (XInvalid) | 0 | 0 |
| 2 | 0 | `gamma_ur` | panic (XInvalid) | 1 | 1 |
| every grid a (0.1 to 1000) | 0 | `gamma_lr` | panic (XInvalid) | 0 | 0 |

Second opinion: SciPy (`scipy.special.gammainc`) returns exactly 0 at all 12
points. The value needs no computation: P(a, x) = gamma(a, x) / Gamma(a) and
gamma(a, 0) is the integral of t^(a-1) e^-t from 0 to 0 (DLMF 8.2.1), so
P(a, 0) = 0 and Q(a, 0) = 1 - P(a, 0) = 1 for every a > 0.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; at x = 0 `mpmath.gammainc(a, 0, 0,
   regularized=True)` is exactly 0, the identity above.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 1,000 shapes a = 10^u with u uniform in [-2, 9] (a fixed
   linear congruential sequence, seed 20261020; `examples/scan.rs` in the
   reproduction project): `gamma_lr(a, 0)` and `gamma_ur(a, 0)` panic for
   all 1,000; `gamma_lr(a, 5e-324)` returns a value for all 1,000.
5. **Root cause** (from the source, `src/function/gamma.rs`):
   `checked_gamma_lr` and `checked_gamma_ur` return `Err(XInvalid)` when
   `x <= 0.0 || x == f64::INFINITY`, although the functions are defined (0
   and 1, or 1 and 0) at both ends; `gamma_lr` and `gamma_ur` are
   `checked_...(a, x).unwrap()`. Accepting x = 0 (and +inf) with those
   values would fix both, as statrs's own distributions already do before
   calling them.

## Notes

The error variant's message ("x is infinite, zero or less than zero")
shows the rejection is deliberate, but the regularized functions are
well-defined there, and every caller that integrates from 0 meets x = 0.
