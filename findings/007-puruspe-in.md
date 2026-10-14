# 007: puruspe::In returns inf or an arbitrary number where I_n(x) is tiny, inf for x > 709.78, and 0.1% errors for x much larger than n

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silent wrong answers: `inf` or `1.1e62` where the true
  value is below 1e-308, `inf` where it is 3.3e306, a relative error of 1.2e-3
  at I_5(700))
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 48 of 200 points of `reference/data/iv.jsonl` (25 of
  120, tolerance 1e-13) and `reference/data/iv_large.jsonl` (23 of 80,
  tolerance 2e-13); the exact lists are `IN` and `IN_LARGE` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Three separate
  bugs, all for n >= 2 except the second:
  1. **Tiny results** (26 points: all 25 of `iv`, and (500, 50) of
     `iv_large`): wherever the true I_n(x) is below 2.2e-308, `In` returns
     `inf` or an arbitrary finite number (e.g. I_500(1) = 1.1e62, true
     2.5e-1285). Off the grid (random check below) it is wrong at 2,533 of
     the 2,583 such points, and also returns `inf` at 7 points where the
     true value is a small normal number (1.8e-213 to 2.3e-306, e.g.
     I_233(8.49)).
  2. **x > 709.78** (13 points of `iv_large`: every order at x = 710, and
     n = 100, 200, 500 at x = 720): `inf` where the true value is finite,
     e.g. I_0(710) = 3.3e306 and I_500(720) = 1.2e238. This includes n = 0
     and 1. Off the grid: all 17 random points with x > 709.78 and a finite
     true value.
  3. **x much larger than n** (9 points of `iv_large`: n = 2 at x = 100 to
     700, n = 5 at x = 200 to 700, n = 10 at x = 500 and 700): relative
     errors from 2.1e-10 (I_2(100)) to 1.2e-3 (I_5(700)). Off the grid: up
     to 2.3e-4 (I_2(555.3)).

Where none of these applies, `In` is accurate: worst 1.2e-14 on `iv` and
8.3e-15 on `iv_large` among passing points.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::In(500, 1.0)); // 1.0952185795074787e62, should be 0 (true 2.5e-1285)
    println!("{:?}", puruspe::In(200, 0.5)); // inf, should be 0 (true 4.9e-496)
    println!("{:?}", puruspe::In(233, 8.492357210415177)); // inf, should be 2.348621565598479e-306
    println!("{:?}", puruspe::In(0, 710.0)); // inf, should be 3.345334558619656e306
    println!("{:?}", puruspe::In(5, 700.0)); // 1.504375359287016e302, should be 1.5025025321928688e302
    println!("{:?}", puruspe::In(2, 555.3081409425245)); // 2.4793472202396365e239, should be 2.479914424874581e239
}
```

## Expected

| Input (n, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (500, 1) | 1.0952185795074787e62 | 0 (2.505e-1285) | 0 | inf |
| (200, 0.5) | inf | 0 (4.912e-496) | 0 | inf |
| (233, 8.492357210415177) | inf | 2.348621565598479e-306 | 2.3486215655982102e-306 | inf |
| (0, 710) | inf | 3.345334558619656e306 | 3.3453345586196545e306 | inf |
| (500, 720) | inf | 1.230559498009891e238 | 1.2305594980098783e238 | inf |
| (100, 720) | inf | 7.098312758841151e307 | inf | inf |
| (5, 700) | 1.504375359287016e302 | 1.5025025321928688e302 | 1.5025025321928705e302 | 1.2e-3 |
| (2, 555.3081409425245) | 2.4793472202396365e239 | 2.479914424874581e239 | 2.4799144248745815e239 | 2.3e-4 |

Rows 3 and 8 are off the grid (random check). "0 (v)" means the true value v
rounds to 0 as a double. All true values are `mpmath.besseli` at 50 digits.

## Notes

The true values come from `mpmath.besseli`, tested in
`reference/tests/test_functions.py` against the power series (DLMF 10.25.2)
up to I_500(20), the large-x expansion (DLMF 10.40.1) at x = 1000, the
Wronskian with `mpmath.besselk` (DLMF 10.28.2) and, for half-integer orders,
closed forms (DLMF 10.49.8). SciPy agrees with mpmath at every failing grid
point except I_100(720), where SciPy returns `inf` ("both wrong" in the
report); there the power series at 400 digits (scratch script) gives
7.0983127588411507709e307, as mpmath. The same series gives I_500(720) and
I_0(710) as mpmath.

Root cause (from the source, `src/bessel.rs`, `In`, which follows Numerical
Recipes' Miller-recurrence routine for I_n; checked with a scratch copy of the code, not
committed, rule 6):

1. **Tiny results**: for n >= 2, `In` runs Miller's downward recurrence and
   rescales the running values by 2^-512 with puruspe's `ldexp` whenever they
   pass 2^512, including the stored value `ans` that becomes I_n. Once `ans`
   is pushed below the double range, `ldexp` returns **infinity** instead of 0
   (the same `ldexp` bug as [finding 006](006-puruspe-jn.md)); further
   rescalings turn that infinity into arbitrary finite numbers. With a correct
   `ldexp` (`x * 2f64.powi(e)`) every one of the 26 grid points gives 0, but
   the 7 off-grid points with a small normal result then give 0 instead of
   their true value: `ans` underflows before the final normalisation
   `ans * i0(x) / bi` brings it back into range. A full fix keeps the
   exponent of `ans` separately (or normalises before it underflows).
2. **x > 709.78**: `i0` and `i1` compute `exp(x) * P(z) / (sqrt(x) Q(z))` for
   x >= 15; `exp(x)` overflows to `inf` above x = 709.78 although the
   quotient is finite up to x = 713.99 (I_0) and beyond for larger n.
   `In` for n >= 2 multiplies by `i0(x)` and so inherits the `inf`
   (scratch: `i0(709) = 1.23e306`, `i0(710) = inf`). A fix would form
   `exp(x / 2) * (P / (sqrt(x) Q)) * exp(x / 2)`.
3. **x much larger than n**: the recurrence starts at order
   2 (n + floor(sqrt(200 n))), which is too low when x is much larger than n
   (Miller's method needs a start well above x for I). In the scratch copy,
   starting 2x + 50 orders higher leaves only the 13 points of bug 2 on the
   grid (worst remaining error 9.3e-15), so the start order alone explains
   the 9 accuracy failures.

Random check (scratch script, not committed, rule 6): 3,000 points with n
uniform on 2..1000 and x = 10^U(-10, 3) (Python `random.seed(20261014)`,
drawn before the 3,000 `Inu_Knu` points of finding 008), each compared with
`mpmath.besseli` at 30 digits under the dataset tolerances (1e-13 for
x <= 20, 2e-13 above). Results: of 2,583 points with a true value below
2.2e-308, `In` returns `inf` at 661 and a wrong finite number at 1,872 (0 at
the other 50). Of 401 points with a normal true value, 24 give `inf` (17 with
x > 709.78, 7 small results listed above) and 2 are inaccurate (I_2(555.3)
2.3e-4, I_8(598.4) 2.4e-8); all 16 points whose true value overflows give
`inf` (accepted).

## Approval checklist

1. Reproduced: the program above was run in a fresh `cargo new` project with
   `puruspe = "=0.4.4"` and printed the values in its comments.
2. Arbiter: mpmath at 50 and 100 digits agree at every point (the generator
   checks); `mpmath.besseli` is tested against the power series, the large-x
   expansion, the Wronskian and closed forms.
3. Second opinion: SciPy agrees with mpmath at 47 of the 48 failing grid
   points; at I_100(720) the power series at 400 digits does.
4. Off-grid: the random check above.
5. Root cause: all three from the source, each confirmed in a scratch copy
   (correct `ldexp`; `i0(710) = inf`; higher start order).
