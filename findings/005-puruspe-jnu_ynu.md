# 005: puruspe::Jnu_Ynu panics for x >= nu + 9845, and loses accuracy as x grows

- **Crate / version:** puruspe 0.4.4
- **Severity:** medium (panic instead of a value; the docs warn "Also panics if
  `x` is too large"), plus a low-severity accuracy loss
- **Status:** approved by owner (2026-10-14)
- **Affected inputs:** 29 of 80 points of `reference/data/yv_large.jsonl`
  (tolerance 5e-10); the exact list is `JNU_YNU_LARGE` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Only the Y_nu
  result is checked here; `besseljy`, `CachedJnuYnu` and the J_nu result share
  the code path below.
  1. **Panic** (28 points): every grid point with x >= 1e4, except
     (200, 1e4) and (500, 1e4), panics with "x too large in besseljy; try
     asymptotic expansion". Off the grid the boundary is x = nu + 9845
     (checked at nu = 0, 500, 1000: x = nu + 9844 returns a value, nu + 9845
     panics).
  2. **Accuracy loss** (1 grid point): Y_500(1e4) has relative error 1.4e-9,
     7.8e2 times the unavoidable loss there (condition number 1.7e4 times
     2^-53). Off the grid, a random check (below) found errors up to 4e-7
     (Y_100(8250), near a zero of Y); relative to the unavoidable loss the
     error grows with x, to about 1.7e3 times it near x = 1e4.

`puruspe::Yn` (integer order) is accurate on the same data (worst 5.5e-11 at
x = 1e6, within the unavoidable loss), so integer orders have a workaround.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::Jnu_Ynu(100.0, 8250.0).1); // -1.4526611486089671e-5, should be -1.4526617240227674e-5
    println!("{:?}", puruspe::Jnu_Ynu(500.0, 1e4).1); // 4.095131064348857e-3, should be 4.095131070250287e-3
    println!("{:?}", puruspe::Jnu_Ynu(0.0, 9845.0).1); // panics: x too large in besseljy
}
```

## Expected

| Input (nu, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (0, 1e4) | panic | 3.6478055589866058e-3 | 3.647805558986606e-3 | inf |
| (1, 1e6) | panic | -3.3104337672417626e-4 | -3.3104337672417637e-4 | inf |
| (0, 9845) | panic | (not on the grid) | | inf |
| (500, 1e4) | 4.095131064348857e-3 | 4.095131070250287e-3 | 4.095131070257228e-3 | 1.4e-9 |
| (100, 8250) | -1.4526611486089671e-5 | -1.4526617240227674e-5 | -1.4526617239389486e-5 | 4.0e-7 |
| (1, 7250) | 9.740269244463216e-5 | 9.740268412194151e-5 | 9.740268412194147e-5 | 8.5e-8 |

The last two rows are off the grid (random check); their mpmath values were
computed at 50 digits with the same function as the dataset.

## Notes

The true values come from `mpmath.bessely`, tested in
`reference/tests/test_functions.py` against the Wronskian
J_{n+1} Y_n - J_n Y_{n+1} = 2/(pi x) (DLMF 10.5.2) at x up to 1e6, the upward
recurrence (DLMF 10.6.1) up to Y_500(20), and the Hankel expansion
(DLMF 10.17.4) at x = 1e6.

Root cause (from the source, `src/bessel.rs`, `besseljy`): the function is
Steed's method as in Numerical Recipes (`bessjy`). Its first continued
fraction (CF1, for J_nu'/J_nu) needs about x - nu + 150 iterations when
x > nu, and the loop stops at `MAXIT = 10000` with
`panic!("x too large in besseljy; try asymptotic expansion")`. A Python
emulation of that loop gives the boundary x = nu + 9845, which matches the
Rust results. Numerical Recipes intends the method for moderate x and
suggests the asymptotic expansion beyond; puruspe has no such branch. The
accuracy loss fits the same loop: rounding errors accumulate over thousands
of CF1 steps. A fix would use the Hankel expansion (DLMF 10.17.3/10.17.4) for
large x, as Boost.Math does (`bessel_jy`, asymptotic branch); not attempted
here.

Random check (scratch script, not committed, rule 6): 3,000 points with
x = 10^U(0, 4) and nu drawn from the grid orders or uniform on [0, 600]
(Python `random.seed(20261014)`), plus x = 5000..20000 step 250 for nine
orders; every result compared with `mpmath.bessely` at 30 digits.
