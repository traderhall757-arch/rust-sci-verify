# 034: puruspe gammq returns 0 in the upper tail and 1 for tiny x, gammp 0 for tiny x, and gammq loses digits through ln_gamma

- **Crate / version:** puruspe 0.4.4 (`puruspe::gammq`, the regularized upper
  incomplete gamma function Q(a, x); `puruspe::gammp`, the lower one P(a, x),
  for its small-x shortcut)
- **Severity:** high (silently wrong answers: `gammq(1, 50)` is 0 instead of
  e^-50 = 1.9e-22, `gammq(0.01, 1e-20)` is 1 instead of 0.365)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 180 of the 538 points of
  `reference/data/gamma_sf.jsonl` for `gammq` (Q(a, x) with rate 1,
  tolerance 1e-12), and 20 of the 538 points of `gamma_cdf` for `gammp`
  (below the groups). The exact lists are `GAMMQ` and `GAMMP_CDF_TINY_X` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Four groups for
  `gammq`, by cause (checklist step 5):
  1. **Upper tail returned as 0** (63 points, every shape on the grid,
     0.01 to 1e6): every x > a + 20 (sqrt(a) + 1) gives exactly 0. The true
     values there go up to 1.6e-21 on the grid (1.1e-15 off it) and are
     normal doubles at all but one of these points.
  2. **Tiny x returned as 1** (6 points, a = 0.01, 0.1, 0.5): every
     x < a 2^-52 gives exactly 1, but for small a, 1 - Q(a, x) is about
     x^a / Gamma(a + 1), which is not small: Q(0.01, 1e-20) = 0.365
     (174% off), Q(0.1, 1e-20) = 0.989 (1.1% off).
  3. **Digits lost below a + 1** (41 points, a = 0.01 to 50): Q is formed
     as 1 - P from the series, and P carries `ln_gamma`'s error
     ([033](033-puruspe-ln-gamma.md)); the subtraction multiplies its
     relative error by P / Q. Up to 8.1e-9 (a = 0.01, x = 1, Q = 0.0022).
  4. **Digits lost above a + 1** (70 points, a = 0.01 to 50): the
     continued fraction's prefactor uses the same `ln_gamma`, so its
     absolute error becomes a relative error in Q: 5.2e-12 at a = 1
     (`gammq(1, 10)` is not e^-10 because `ln_gamma(1)` = -5.2e-12), up to
     4.9e-11 at a = 50.

The other 43 failures on `gamma_sf` (a >= 100, 0.2a <= x <= a + 20
(sqrt(a) + 1)) come from the quadrature `gammpapprox` that `gammq` shares
with `gammp`: 23 NaN and 20 accuracy losses up to 1.2e-9. They are listed
under [finding 002](002-puruspe-gammp.md) (`GAMMQ_APPROX`).

`gammp` has the same two shortcuts. Its upper one (P = 1 for
x > a + 20 (sqrt(a) + 1)) is harmless, since 1 - P is below 2^-53 there.
Its lower one is not: it returns 0 for every x < a 2^-52, where
P(a, x) is about x^a e^-x / Gamma(a + 1), so the answer is wrong unless the
true value underflows. 20 of the 538 points of
`reference/data/gamma_cdf.jsonl` (list `GAMMP_CDF_TINY_X`, a = 0.01 to 10,
x = 1e-300 to 1e-20) fail this way: `gammp(0.01, 1e-20)` = 0 instead of
0.635, `gammp(1, 1e-20)` = 0 instead of 1e-20. The `gammainc` grid starts
at x = 0.001 and never reached it. The other 155 `gamma_cdf` failures of
`gammp` are finding 002's.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::gammq(1.0, 50.0)); // should be 1.9287498479639178e-22 (= e^-50)
    println!("{:?}", puruspe::gammq(10.0, 100.0)); // should be 1.1253473960842733e-31
    println!("{:?}", puruspe::gammq(1e4, 13000.0)); // should be 4.712347177485643e-166
    println!("{:?}", puruspe::gammq(0.01, 1e-20)); // should be 0.36544207945100327
    println!("{:?}", puruspe::gammq(0.1, 1e-20)); // should be 0.9894886299388822
    println!("{:?}", puruspe::gammq(0.01, 1.0)); // should be 0.00221623462322799
    println!("{:?}", puruspe::gammq(1.0, 10.0)); // should be 4.5399929762484854e-05 (= e^-10)
    println!("{:?}", puruspe::gammq(50.0, 57.071067811865476)); // should be 0.15792211722597285
    println!("{:?}", puruspe::gammq(100.0, 50.0)); // should be 0.9999999996799934 (finding 002)
    println!("{:?}", puruspe::gammq(1e6, 1001000.0)); // should be 0.15865521363165971 (finding 002)
    println!("{:?}", puruspe::gammp(0.01, 1e-20)); // should be 0.6345579205489967
    println!("{:?}", puruspe::gammp(1.0, 1e-20)); // should be 1e-20
}
```

Output (release and debug mode identical): `0.0`, `0.0`, `0.0`, `1.0`,
`1.0`, `0.0022162346053057336`, `4.539992976272241e-5`,
`0.15792211721819466`, `NaN`, `0.15865521344059652`, `0.0`, `0.0`.

## Expected

| (a, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (1, 50) | 0.0 | 1.9287498479639178e-22 | 1.9287498479639183e-22 | 1 |
| (10, 100) | 0.0 | 1.1253473960842733e-31 | 1.1253473960842808e-31 | 1 |
| (1e4, 13000) | 0.0 | 4.712347177485643e-166 | 4.7123471774867427e-166 | 1 |
| (0.01, 1e-20) | 1.0 | 0.36544207945100327 | 0.36544207945100327 | 1.7 |
| (0.1, 1e-20) | 1.0 | 0.9894886299388822 | 0.9894886299388822 | 1.1e-2 |
| (0.01, 1) | 0.0022162346053057336 | 0.00221623462322799 | 0.00221623462322799 | 8.1e-9 |
| (1, 10) | 4.539992976272241e-5 | 4.5399929762484854e-5 | 4.539992976248486e-5 | 5.2e-12 |
| (50, 57.071067811865476) | 0.15792211721819466 | 0.15792211722597285 | 0.15792211722597285 | 4.9e-11 |
| `gammp` (0.01, 1e-20) | 0.0 | 0.6345579205489967 | 0.6345579205489967 | 1 |

Second opinion: SciPy (`scipy.special.gammainc`) agrees with mpmath within
5.1e-14 at the 20 `gammp` points. SciPy (`scipy.special.gammaincc`) agrees
with mpmath within 2.6e-14 at all 180 `gammq` points but one, (5, 750), where the true value
2.5e-316 is subnormal and SciPy flushes it to 0 as well (verdict "both
wrong"). Q(1, x) = e^-x is a closed form for the rows with a = 1.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `puruspe = "=0.4.4"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double; the `gamma_sf` route (series DLMF 8.7.1,
   continued fraction DLMF 8.9.2) is tested against `mpmath.gammainc`, the
   finite sum DLMF 8.4.10 and closed forms in
   `reference/tests/test_functions.py` (finding 012).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 2,000 random points (Python `random.Random(20261030)`):
   a = 10^u with u uniform in [-2, 6]; x = a 10^v with v uniform in
   [-25, -14] (15%) or [-3, 0.5] (30%), x = a + z sqrt(a) with z uniform in
   [-10, 45] (40%), or x uniform in [0, 745] (15%); 6 significant digits;
   reference by the dataset's route at 50 and 100 digits, which agree at
   every point. Against the tolerance 1e-12:
   - x > a + 20 (sqrt(a) + 1): 401 of 449 fail, each returning 0 for a true
     value from 2.2e-316 to 1.1e-15; the 48 that pass have true values
     below the smallest normal double;
   - x < a 2^-52: 51 of 300 fail (returning 1 for true values from 0.454
     to 1 - 1e-16), every failing shape is at most 0.664;
   - series region (x < a + 1, a < 100, or x < 0.2a for a >= 100): 177 of
     684 fail, worst 7.1e-9 (a = 0.0106268, x = 0.962733);
   - continued-fraction region (a < 100, x >= a + 1, below the cut-off):
     242 of 243 fail, worst 9.3e-11;
   - quadrature region (a >= 100, a - 1 < x below the cut-off): 175 of 175
     fail, worst 1.8e-9, and NaN at all 149 points with
     0.2a <= x <= a - 1 (finding 002);
   - `gammp` with x < a 2^-52: 500 random points (`random.Random(20261031)`,
     a = 10^u with u uniform in [-2, 3], x = 10^w with w uniform from -300
     to log10(a 2^-52), 6 significant digits; same reference): at all 227
     points where the true P is a normal double, `gammp` returns 0 (true
     values up to 0.32; the largest failing shape is 9.41); at the 273
     where P is subnormal or 0, 6 fail (flushed subnormals).
5. **Root cause** (from the source, `src/gamma.rs`), confirmed with a copy of
   `gammq`, `gser`, `gcf` and `gammpapprox` (a scratch crate, not committed)
   whose switches turn each cause off; with all switches off it matches
   `puruspe::gammq` bit for bit at all 538 grid points.
   - **Group 1**: `if x > a + 20.0 * (a.sqrt() + 1.0) { return 0f64; }`.
     This is a shortcut for P = 1 (Q is below about 1e-15 there), copied
     into `gammq`, where it discards the whole answer.
   - **Group 2**: `if x < a * EPS { return 1f64; }`, the mirror shortcut
     for P = 0, which is wrong for small a because P(a, x) is about
     x^a / Gamma(a + 1).
   - **Groups 3 and 4**: `gser` and `gcf` compute their prefactor
     exp(-x + a ln x - ln Gamma(a)) with puruspe's 7-term `ln_gamma`,
     whose absolute error up to 1.3e-10 (finding 033) is a relative error
     in the result; in group 3, 1 - P then amplifies it.
   - `gammp`: `if x < a * EPS { return 0f64; }`, the same shortcut with
     the same flaw. The copy (which also matches `puruspe::gammp` bit for
     bit on `gamma_cdf` and `gammainc`) without it passes all 20 grid
     points and all 227 random ones above.
   - With both cut-offs removed, `ln_gamma` replaced by libm's `lgamma`
     (about 1 ulp) and the -700 shortcut below dropped, the copy fails 40
     grid points: the 23 NaN and 17 in
     the quadrature region (finding 002), none in groups 1 to 4. Off the
     grid (same 2,000 points) the copy fails 0 points in the series,
     continued-fraction and small-x regions; with x above the cut-off it
     fails 0 for a < 100 and 128 of 231 for a >= 100, where the
     quadrature takes over (up to 3.7e-5 far in the tail; finding 002's
     lead). A third, smaller cut-off (`log_result < -700.0` returns 0 in
     `gser` and `gcf`) flushes results below about 1e-304 to 0; on the grid
     it only matters at points group 1 already covers.

## Notes

Leads for a fix: drop both shortcuts from `gammq` (or replace them with
the exact limits: Q = 1 - x^a / Gamma(a + 1) for tiny x, and the continued
fraction in the upper tail, which converges fast there), compute the
prefactor with a better `ln_gamma` or the Lanczos/Temme form that Boost uses
(`regularised_gamma_prefix`), and take Q directly from the continued
fraction or Temme's expansion instead of 1 - P where P is close to 1.
`gammq` is not used inside puruspe; `invgammp` uses `gammp`.
