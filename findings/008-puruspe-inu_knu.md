# 008: puruspe::Inu_Knu gives NaN for I where I_nu(x) underflows, inf for x above about 714, and loses accuracy

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silent wrong answer: `inf` where the true I_nu(x) is
  1.2e238), plus NaN instead of 0 (medium) and an accuracy loss up to 3.7e-7
  (low)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** the I_nu result (the first element of the tuple) at 49
  of 296 points: 25 of 120 in `reference/data/iv.jsonl` (tolerance 1e-13),
  9 of 80 in `reference/data/iv_large.jsonl` (2e-13) and 15 of 96 in
  `reference/data/iv_nonint.jsonl` (2e-13); the exact lists are `INU_IV`,
  `INU_IV_LARGE` and `INU_IV_NONINT` in `crates/check-puruspe/src/expected.rs`,
  the report is [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md).
  `besselik` and the cached variants share the code path. The K_nu result is
  not checked here.
  1. **NaN where I underflows** (23 points: 19 of `iv`, 4 of `iv_nonint`):
     e.g. I_100(1e-5), I_500(20), I_200.75(0.1); the right answer is 0. Off
     the grid: NaN at 1,317 of the 1,583 random points whose true I is below
     2.2e-308, and never where it is normal.
  2. **inf for x above about 714** (3 points of `iv_large`: nu = 100, 200,
     500 at x = 720): the true values 7.1e307, 7.4e298 and 1.2e238 are
     finite. Off the grid: all 6 random points with a finite true value and
     x > 714 (x from 716.6 to 804.8).
  3. **Accuracy loss** (23 points), low severity:
     - large orders (nu >= 100): 6 points of `iv` (nu = 100, 200 at
       x = 5, 10, 20; 3.9e-13 to 8.8e-13), 6 of `iv_large` (nu = 200 at
       x = 50, nu = 500 at x = 100 to 710; up to 1.2e-11 at I_500(100)) and
       7 of `iv_nonint` (nu = 100.5, 200.75 at x = 5 to 50; up to 8.2e-13);
     - small x with a negative fractional part (below): 4 points of
       `iv_nonint` at x = 1e-5 (nu = 0.5, 1.5, 2.7, 10.5; up to 2.0e-10 at
       I_10.5(1e-5)).
     Off the grid: 323 of 1,384 random points with a normal true value are
     outside the tolerance; worst 3.7e-7 (I_5.5048(2.17e-9), small x) and
     1.2e-11 (I_519.25(109.2), large order).

Elsewhere `Inu_Knu`'s I result is within the tolerances (worst passing
4.4e-14 on `iv`, 1.5e-13 on `iv_large`, 1.7e-13 on `iv_nonint`).
`puruspe::In` (integer order) has its own bugs: finding 007.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    println!("{:?}", puruspe::Inu_Knu(100.0, 1e-5).0); // NaN, should be 0 (true 8.5e-689)
    println!("{:?}", puruspe::Inu_Knu(500.0, 720.0).0); // inf, should be 1.230559498009891e238
    println!("{:?}", puruspe::Inu_Knu(500.0, 100.0).0); // 3.5900514314256536e-283, should be 3.5900514314696126e-283
    println!("{:?}", puruspe::Inu_Knu(10.5, 1e-5).0); // 1.8350974827260546e-63, should be 1.8350974823518723e-63
    println!("{:?}", puruspe::Inu_Knu(5.504755644692505, 2.173580493901863e-9).0); // 1.5604487831489151e-52, should be 1.5604482105762322e-52
}
```

## Expected

| Input (nu, x) | Crate I | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (100, 1e-5) | NaN | 0 (8.453e-689) | 0 | inf |
| (500, 720) | inf | 1.230559498009891e238 | 1.2305594980098783e238 | inf |
| (100, 720) | inf | 7.098312758841151e307 | inf | inf |
| (500, 100) | 3.5900514314256536e-283 | 3.5900514314696126e-283 | 3.590051431470047e-283 | 1.2e-11 |
| (10.5, 1e-5) | 1.8350974827260546e-63 | 1.8350974823518723e-63 | 1.8350974823518745e-63 | 2.0e-10 |
| (5.504755644692505, 2.173580493901863e-9) | 1.5604487831489151e-52 | 1.5604482105762322e-52 | 1.5604482105762329e-52 | 3.7e-7 |

The last row is off the grid (random check). All true values are
`mpmath.besseli` at 50 digits.

## Notes

The true values come from `mpmath.besseli`, tested as described in
[finding 007](007-puruspe-in.md) (power series, large-x expansion, Wronskian
with `besselk`, half-integer closed forms). SciPy agrees with mpmath at every
failing grid point except I_100(720), where it returns `inf` ("both wrong" in
the report) and the power series at 400 digits agrees with mpmath.

Root cause (from the source, `src/bessel.rs`, `besselik`, which is Temme's
method as in Numerical Recipes `besselik`; each item checked with a scratch
copy of the code that prints its intermediate values, not committed, rule 6):

1. **NaN**: an unnormalised I is recurred downwards from `ril = FPMIN`
   (about 1e-292) over `nl = round(nu)` orders. When I_xmu(x) / I_nu(x)
   exceeds about 1e600, `ril` overflows to `inf` and `f = ripl / ril` is NaN
   (scratch at (500, 20): `ril = inf`, `f = NaN`). The same overflow as the J
   result of [finding 005](005-puruspe-jnu_ynu.md).
2. **inf**: I_xmu is found from the Wronskian,
   `rimu = (1/x) / (f K_xmu - K_xmu')`. For x above about 708, K_xmu(x)
   (about e^-x) is subnormal, and the true I_xmu(x) is above the double range
   from about x = 714 even when I_nu(x) for the requested order is not.
   Scratch at (500, 720): `K_xmu = 9.5e-315`, denominator `1.9e-314`,
   `rimu = inf`. A fix would scale K and I by e^x and e^-x.
3. **Accuracy loss, large orders**: the recurrence factor `fact` is formed
   as nu / x minus `xi` (1 / x) repeated `nl` times instead of as
   (nu - k) / x, so its rounding errors accumulate. Confirmed: the scratch
   copy with `fact = (nu - k) * xi` passes all 19 large-order points (worst
   remaining error 1.0e-14 on `iv_large`, 9.6e-15 on `iv`). The J result of
   finding 005 builds `fact` the same way (its lead for the J accuracy loss).
4. **Accuracy loss, small x**: `nl` rounds nu to the nearest integer, so the
   fractional order xmu is in [-1/2, 1/2). For negative xmu, I_xmu(x) is
   nearly proportional to K_xmu(x) at small x (it contains the term
   (2/pi) sin(|xmu| pi) K_|xmu|(x)), and the Wronskian denominator
   `f K - K'` cancels: scratch at (0.5, 1e-5): `f K = -1.98164e7`,
   `K' = -1.98168e7`, difference 396, a loss of about 5e4 x 2^-53 (observed
   error 2.9e-12). The loss grows like x^(-2|xmu|). With the `fact` fix
   above, (2.7, 1e-5) passes and (0.5, 1.5, 10.5 at 1e-5) still fail
   (2.9e-12 to 6.5e-12). The x < 2 branch of the J result has the same
   structure, which probably explains finding 005's small-x loss for
   negative xmu (not checked).

Random check (scratch script, not committed, rule 6): 3,000 points with nu
drawn from the grid orders or uniform on [0, 600] (half each) and
x = 10^U(-10, 3) (Python `random.seed(20261014)`, drawn after the 3,000 `In`
points of finding 007), compared with `mpmath.besseli` at 30 digits under
the dataset tolerances (1e-13 for integer nu and x <= 20, 2e-13 otherwise).
Results: of 1,583 points whose true I is below 2.2e-308, 1,317 give NaN,
1 a wrong subnormal (6.3e-12 relative to 2.2e-308), the rest 0. Of 1,384
points with a normal true value, 6 give `inf` (x > 714), 323 are inaccurate
(above); all 33 points whose true value overflows give `inf` (accepted). No
panics.

## Approval checklist

1. Reproduced: the program above was run in a fresh `cargo new` project with
   `puruspe = "=0.4.4"` and printed the values in its comments.
2. Arbiter: mpmath at 50 and 100 digits agree at every point (the generator
   checks); `mpmath.besseli` is tested against independent routes (finding
   007).
3. Second opinion: SciPy agrees with mpmath at 48 of the 49 failing grid
   points; at I_100(720) the power series at 400 digits does.
4. Off-grid: the random check above.
5. Root cause: all four from the source, confirmed with scratch copies
   (intermediate values; the `fact` change).
