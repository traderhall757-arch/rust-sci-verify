# 006: puruspe::Jn returns inf or an arbitrary number where J_n(x) underflows

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silent wrong answer: `inf`, `1.8e90` or `4.6e-3` where
  the true value is below 1e-308)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 27 of 200 points of `reference/data/jv.jsonl` (26 of
  120, tolerance 1e-13) and `reference/data/jv_large.jsonl` (1 of 80,
  tolerance 5e-10); the exact lists are `JN` and `JN_LARGE` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Every failing
  point is one where the true J_n(x) is below the smallest normal double
  (2.2e-308), so the right answer is 0 (or a subnormal). Off the grid
  (random check below) `Jn` is wrong at 2,438 of 2,491 such points and at
  none of the 509 points where the true value is a normal double. In
  practice: n >= 2 and x well below n, e.g. J_500(x) for x <= 50, J_200(x)
  for x <= 2, J_1000(100).

Where J_n(x) is a normal double, `Jn` is accurate on both datasets (worst
1.2e-14 on `jv`, 6.7e-11 on `jv_large`, both within the unavoidable loss).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::Jn(500, 1.0)); // 0.004594918722316429, should be 0 (true value about 1e-1300)
    println!("{:?}", puruspe::Jn(500, 20.0)); // 1.7708659106201395e90, should be 0
    println!("{:?}", puruspe::Jn(200, 0.5)); // inf, should be 0
    println!("{:?}", puruspe::Jn(1000, 100.0)); // 1.826802375194446e-188, should be 0 (true 1.9e-870)
}
```

## Expected

| Input (n, x) | Crate | mpmath (true value) | SciPy |
|---|---|---|---|
| (500, 1) | 4.594918722316429e-3 | 0 (about 1e-1300) | 0 |
| (500, 20) | 1.7708659106201395e90 | 0 | 0 |
| (200, 0.5) | inf | 0 | 0 |
| (500, 50) | inf | 0 (2.19e-436) | 0 |
| (1000, 100) | 1.826802375194446e-188 | 0 (1.90e-870) | 0 |
| (300, 30) | 1.0388021531643552e-262 | 1.03880215316435e-262 | 1.0388021531643144e-262 |

The last two rows are off the grid; the last one shows a small but normal
J_n(x) in the same branch, which `Jn` gets right (relative error 5e-15).
"0 (v)" means the true value v rounds to 0 as a double; all such values
were computed with `mpmath.besselj` at 30 digits.

## Notes

The true values come from `mpmath.besselj`, tested in
`reference/tests/test_functions.py` against the power series (DLMF 10.2.2)
up to J_500(20) and at J_500(1e-10), the Hankel expansion (DLMF 10.17.3) at
x = 1e6, and, through the Wronskian test, against `mpmath.bessely`. SciPy
agrees with mpmath at every failing grid point (it returns 0).

Root cause (from the source): for x < n, `Jn` (`src/bessel.rs`) uses
Miller's downward recurrence as in Numerical Recipes (`bessj`). To avoid
overflow it rescales the running values by 2^-512 with puruspe's own
`ldexp` (`src/utils.rs`) whenever they grow past 2^512. Once the value that
will become J_n has been stored (`ans`), further rescaling pushes it below
the double range, and `ldexp` then returns **infinity** instead of 0: its
range check `if !(0..=0x7ff).contains(&new_exponent)` returns `+-inf` for
both overflow and underflow, although its comment says "infinity or zero".
A later rescaling of that infinity gives a finite number (exponent field
0x7ff minus 512 is a valid exponent), so the result is `inf` or an
arbitrary finite value depending on how many rescalings follow. The same
`ldexp` also mishandles a new exponent field of exactly 0 or 0x7ff, and
subnormal inputs.

Confirmed in a scratch program (not committed, rule 6): the same `Jn` code
with only `ldexp(x, e)` replaced by `x * 2f64.powi(e)` returns 0 at every
underflow point above and is bit-identical to puruspe where J_n(x) is
normal (checked at (100, 5), (200, 20), (500, 200)). A fix would be that
change, or a correct `ldexp` (return a signed zero or a subnormal on
underflow). `In` (modified Bessel, integer order) uses the same rescaling
with the same `ldexp`; it is checked in the next dataset (`iv`).

Random check (scratch script, not committed, rule 6): 3,000 points with
n uniform on 2..1000 and x = 10^U(-10, 3.5) (Python `random.seed(20261014)`,
drawn before the 3,000 `Jnu_Ynu` points of finding 005), each compared with
`mpmath.besselj` at 30 digits under the dataset tolerances (1e-13 for
x <= 20, 5e-10 above). Results: 2,491 points have a true value below
2.2e-308; `Jn` returns `inf` at 674 and a wrong finite number at 1,764, and
0 at the other 53. No failure with x > n (the upward-recurrence branch) and
none where the true value is normal; the largest x / n among failures is
0.27.

## Approval checklist

1. Reproduced: the program above was run in a fresh `cargo new` project
   with `puruspe = "=0.4.4"`.
2. Arbiter: mpmath at 50 and 100 digits agree at every point (the generator
   checks); `mpmath.besselj` is tested against the power series and the
   Hankel expansion.
3. Second opinion: SciPy agrees with mpmath at all 27 failing grid points.
4. Off-grid: the random check above.
5. Root cause: the `ldexp` underflow branch, confirmed by the scratch fix.
