# 005: puruspe::Jnu_Ynu panics for x above about nu + 9832, gives NaN for J where J underflows, and loses accuracy

- **Crate / version:** puruspe 0.4.4
- **Severity:** medium (panic instead of a value; the docs warn "Also panics if
  `x` is too large"; NaN instead of 0), plus a low-severity accuracy loss
- **Status:** approved by owner (2026-10-14) for the Y_nu result; the J_nu
  part (section "The J_nu result" below) approved (Claude, checklist
  complete, 2026-10-14)
- **Affected inputs:** 29 of 80 points of `reference/data/yv_large.jsonl`
  (tolerance 5e-10); the exact list is `JNU_YNU_LARGE` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). This list
  checks the Y_nu result; the J_nu result is in "The J_nu result" below.
  `besseljy` and `CachedJnuYnu` share the code path.
  1. **Panic** (28 points): every grid point with x >= 1e4, except
     (200, 1e4) and (500, 1e4), panics with "x too large in besseljy; try
     asymptotic expansion". Off the grid the first panic is at about
     x = nu + 9832 and every x above about nu + 9851 panics; in between,
     panics and values alternate (scan of x = nu + 9800 to nu + 9900 in steps
     of 0.05 at nu = 0, 500, 1000: first panic at nu + 9832.4, 9833.95,
     9832.35; last value at nu + 9850.85, 9849.1, 9845.9). Corrected
     2026-10-14 with the J_nu part: the first version said "x = nu + 9845
     exactly", from checks at nu + 9844 and nu + 9845 only.
  2. **Accuracy loss** (1 grid point): Y_500(1e4) has relative error 1.4e-9,
     7.8e2 times the unavoidable loss there (condition number 1.7e4 times
     2^-53). Off the grid, a random check (below) found errors up to 4e-7
     (Y_100(8250), near a zero of Y); relative to the unavoidable loss the
     error grows with x, to about 1.7e3 times it near x = 1e4.

`puruspe::Yn` (integer order) is accurate on the same data (worst 5.5e-11 at
x = 1e6, within the unavoidable loss), so integer orders have a workaround.

## The J_nu result

Added with the `jv` and `jv_large` datasets (same grids as `yv` and
`yv_large`, tolerances 1e-13 and 5e-10). 55 of 200 points fail; the lists
are `JNU_JV` (25 points) and `JNU_JV_LARGE` (30) in
`crates/check-puruspe/src/expected.rs`.

1. **Panic** (28 points of `jv_large`): the same points and the same cause
   as for Y_nu (item 1 above), with the same off-grid boundary (one call
   returns both).
2. **NaN where J underflows** (19 points of `jv`): 19 of the 27 grid points
   where the true J_nu(x) is below 2.2e-308 (e.g. J_100(1e-5), J_500(20))
   give NaN; the right answer is 0. The other 8 give 0. (The Y_nu result at the same
   points is -inf or NaN, both accepted at an overflow point.) Off the grid:
   NaN at 1,223 of the 1,464 random points whose true J is below 2.2e-308,
   and never where it is normal.
3. **Accuracy loss**, low severity:
   - 6 points of `jv`, nu = 100 and 200 at x = 5, 10, 20: relative error
     3.8e-13 to 8.9e-13, 4 to 9 times the tolerance and 17 to 40 times the
     unavoidable loss (condition number about nu, times 2^-53). At
     J_200(5) SciPy is also wrong (it returns 0 for 4.76e-296), so the
     report lists that point as "both wrong".
   - 2 points of `jv_large`: J_200(1e4) has error 1.9e-8 (SciPy 2.0e-10),
     J_500(1e4) 5.1e-10, from the same long CF1 loop as Y.
   - Off the grid, non-integer orders at small x (x < 2) are worse: errors
     up to 1.3e-9 (J_39.549389116250055(8.17756317156625e-7)), mostly where
     the fractional part xmu = nu - round(nu) is negative (46 of 61 such
     random points fail, worst 1.3e-9; for xmu >= 0, 104 of 776 fail, worst
     1.8e-12).

`puruspe::Jn` (integer order) is accurate where J is a normal double but
has its own bug where J underflows: finding 006.

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::Jnu_Ynu(100.0, 1e-5).0); // NaN, should be 0
    println!("{:?}", puruspe::Jnu_Ynu(200.0, 10.0).0); // 6.967530155398942e-236, should be 6.967530155393544e-236
    println!("{:?}", puruspe::Jnu_Ynu(200.0, 1e4).0); // -3.634005166175012e-4, should be -3.634005234268351e-4
    println!("{:?}", puruspe::Jnu_Ynu(39.549389116250055, 8.17756317156625e-7).0); // 1.4257577369361821e-300, should be 1.4257577351449877e-300
}
```

| Input (nu, x) | Crate J | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (100, 1e-5) | NaN | 0 | 0 | inf |
| (500, 20) | NaN | 0 | 0 | inf |
| (200, 10) | 6.967530155398942e-236 | 6.967530155393544e-236 | 6.967530155393012e-236 | 7.7e-13 |
| (200, 5) | 4.7600100880065223e-296 | 4.76001008800229e-296 | 0 | 8.9e-13 |
| (200, 1e4) | -3.634005166175012e-4 | -3.634005234268351e-4 | -3.634005234983698e-4 | 1.9e-8 |
| (39.549389116250055, 8.17756317156625e-7) | 1.4257577369361821e-300 | 1.4257577351449877e-300 | 1.4257577351450376e-300 | 1.3e-9 |

The last row is off the grid. Second opinion where SciPy is wrong or close
to the tolerance: the power series (DLMF 10.2.2) at 60 digits gives
J_200(5) = 4.76001008800229013...e-296 and J_200(10), J_100(5) as mpmath;
the Hankel expansion (DLMF 10.17.3) at 60 digits gives J_200(1e4) and
J_500(1e4) as mpmath to 20 digits (scratch script). The same two routes are
tested in `reference/tests/test_functions.py`.

Root cause and leads (`src/bessel.rs`, `besseljy`, Numerical Recipes
`bessjy`):
- NaN: the function recurs an unnormalised J downwards from
  `rjl = FPMIN` (about 1e-292) over `nl` orders, about nu - x steps. When
  J_xmu(x) / J_nu(x) exceeds about 1e600 (that is, when J_nu(x) underflows),
  `rjl` and `rjpl` overflow to infinity and `f = rjpl / rjl` is NaN; NaN then
  reaches both results. Numerical Recipes has no overflow guard there
  either. A fix would rescale during the recurrence, or return 0 for J when
  it is known to underflow.
- Accuracy loss for integer nu at x = 5..20 (lead, not tested): the
  recurrence factor `fact` is formed as nu / x minus `xi` repeated `nl`
  (about nu - x) times instead of as (nu - k) / x, so its rounding errors
  accumulate over the recurrence; nl * nu * 2^-53 / x is 2e-13 to 4e-13 at
  the failing points, the same order as the errors.
- Accuracy loss at small x for negative xmu (lead): the x < 2 branch
  (Temme's series for Y_xmu, then J_xmu from the Wronskian
  `rjmu = w / (rymup - f * rymu)`); the cause is not identified.

Random check (scratch, not committed): 3,000 points, nu drawn from the grid
orders or uniform on [0, 600] (half each), x = 10^U(-10, 4), Python
`random.seed(20261014)` (drawn after the 3,000 `Jn` points of finding 006),
compared with `mpmath.besselj` at 30 digits under the dataset tolerances.
1,481 fail: 1,223 NaN (all where J underflows), 3 panics (smallest
x - nu = 9849.8), 255 finite values outside the tolerance (197 with x <= 20,
58 above).

Approval checklist for the J_nu part: (1) the program above was run in a
fresh `cargo new` project; (2) mpmath at 50 and 100 digits agrees at every
point, and `mpmath.besselj` is tested against the power series and the
Hankel expansion; (3) SciPy agrees with mpmath at the NaN points and within
tolerance at the accuracy points except J_200(5), where the power series is
the second route; (4) the random check above; (5) root cause for the panic
and the NaN, leads for the accuracy loss.

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
`panic!("x too large in besseljy; try asymptotic expansion")`. Its stopping
test `|del - 1| <= EPS` is met or missed by rounding near the end, which is
why panics and values alternate between about nu + 9832 and nu + 9851. Numerical Recipes intends the method for moderate x and
suggests the asymptotic expansion beyond; puruspe has no such branch. The
accuracy loss fits the same loop: rounding errors accumulate over thousands
of CF1 steps. A fix would use the Hankel expansion (DLMF 10.17.3/10.17.4) for
large x, as Boost.Math does (`bessel_jy`, asymptotic branch); not attempted
here.

Random check (scratch script, not committed, rule 6): 3,000 points with
x = 10^U(0, 4) and nu drawn from the grid orders or uniform on [0, 600]
(Python `random.seed(20261014)`), plus x = 5000..20000 step 250 for nine
orders; every result compared with `mpmath.bessely` at 30 digits.
