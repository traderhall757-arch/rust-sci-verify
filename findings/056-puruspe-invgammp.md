# 056: puruspe invgammp is wrong in both tails, finite at p = 1, panics for large shapes, and inherits gammp's errors

- **Crate / version:** puruspe 0.4.4 (`puruspe::invgammp(p, a)`, the x with
  P(a, x) = p: the inverse of the regularized lower incomplete gamma
  function, i.e. the quantile of the gamma distribution with rate 1)
- **Severity:** high (silently wrong answers: `invgammp(0.5, 0.01)` is
  2.0e-28 instead of 4.5e-31, `invgammp(1e-50, 2)` is 1.9e-9 instead of
  1.4e-25; a panic "Bad args in gammp" for ordinary inputs such as
  `invgammp(0.1, 1000)`)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 397 of the 432 points of
  `reference/data/gamma_ppf.jsonl` (tolerance 2e-13); every shape on the
  grid, 0.01 to 1e6, fails somewhere. The exact lists are `INVGAMMP`
  (107 points, groups 1, 3, 4 and 5: invgammp's own code) and
  `INVGAMMP_GAMMP` (290 points, groups 2 and 6: errors of the `gammp` it
  calls) in `crates/check-puruspe/src/expected.rs`; the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md). Six groups, by
  cause (checklist step 5):
  1. **p = 1 gives a finite number** (16 points, every shape): it returns
     max(100, a + 100 sqrt(a)) (`invgammp(1, 1)` = 101), but P(a, x) < 1
     for every finite x, so the quantile is +inf (the point is tagged
     `pole`).
  2. **Panic through gammp's NaN** (61 points, a >= 100, p <= 0.4): the
     Halley step calls `gammp`, which returns NaN for 0.2a <= x <= a - 1
     ([finding 002](002-puruspe-gammp.md)); the NaN becomes the next x, and
     the next `gammp(a, NaN)` fails its argument check and panics with
     "Bad args in gammp".
  3. **Absolute stopping test, a < 1** (19 points, a = 0.01, 0.1, 0.5,
     p up to 0.6, where the quantile is far below 1): the loop stops when the step is below
     max(x 2^-52, 2^-52), and for a quantile far below 1 the first step
     already is, so it stops after one step: 5.1x to 453x too large
     (`invgammp(0.5, 0.01)` = 2.0e-28, true 4.5e-31). The start there
     (1.0e-30) is below a 2^-52, where `gammp` returns 0
     ([finding 034](034-puruspe-gammq.md)), so that one step even goes
     the wrong way; with an exact P the single step still leaves 5%
     (4.24e-31).
  4. **Lower tail for a > 1: a start at 1e-3 and 12 steps** (35 points,
     a = 1.5 to 100, p <= 1e-20): for a > 1 the starting guess is clamped
     to at least 1e-3. When the true quantile is far below that, each
     Halley step shrinks x only by a bounded factor (3 for a = 2, where
     P ~ x^2 / 2; a step past 0 is replaced by x/2), and 12 steps do not
     get there (`invgammp(1e-50, 2)` = 1.9e-9 = 1e-3 / 3^12, true 1.4e-25;
     `invgammp(5e-324, 1.5)` = 4.1e-12, true 3.5e-216). Where the density
     at the start underflows to 0 the step is -inf / inf = NaN and the next
     `gammp` call panics (2 points, a = 100, p = 5e-324 and 1e-300).
  5. **Upper tail through P** (37 points, every shape, p >= 1 - 1e-5): the
     equation solved is P(a, x) - p = 0. Near p = 1, P is a double next to
     1 with an absolute error of at least 2^-53, so x is only known to about
     2^-53 / f(x), f the density; it should solve Q(a, x) = 1 - p (exact
     for p >= 1/2) instead. Up to 8.7% off (`invgammp(1 - 1e-15, 0.01)` =
     24.35, true 26.65; `invgammp(1 - 1e-15, 1)` = 34.559, true
     -ln(1e-15) = 34.540).
  6. **gammp's own errors** (229 points, every shape): the iteration
     converges to the root of puruspe's P, not of the true P. 225 points
     lose 2.3e-13 to 2.3e-9 (`ln_gamma` in the series and continued
     fraction prefactor, [finding 033](033-puruspe-ln-gamma.md), for
     a < 100; the quadrature, [finding 002](002-puruspe-gammp.md), for
     a >= 100; `invgammp(0.5, 1)` = 0.6931471805547 instead of
     ln 2 = 0.6931471805599). 4 points (a = 1, p = 1e-300 to 1e-20) are
     2.59x too large: `gammp` returns 0 for x < a 2^-52
     ([finding 034](034-puruspe-gammq.md)), so the steps see P = 0 and
     stop by group 3's test.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // invgammp(p, a): the x with P(a, x) = p
    println!("{:?}", puruspe::invgammp(1.0, 1.0)); // should be inf (P(1, x) = 1 - e^-x < 1 for every finite x)
    println!("{:?}", puruspe::invgammp(0.5, 0.01)); // should be 4.465535018910355e-31
    println!("{:?}", puruspe::invgammp(1e-50, 2.0)); // should be 1.4142135623730952e-25
    println!("{:?}", puruspe::invgammp(0.999999999999999, 1.0)); // should be 34.53957599234088 (= -ln(1 - p))
    println!("{:?}", puruspe::invgammp(1e-20, 1.0)); // should be 1e-20 (to 17 digits)
    println!("{:?}", puruspe::invgammp(0.5, 1.0)); // should be 0.6931471805599453 (= ln 2)
    println!("{:?}", puruspe::invgammp(0.99, 0.01)); // should be 0.26505255025158186
    println!("{:?}", puruspe::invgammp(0.5, 100.0)); // should be 99.66686491931549
    let r = std::panic::catch_unwind(|| puruspe::invgammp(0.1, 1000.0)); // should be 959.6939327288333
    println!("{r:?}");
    let r = std::panic::catch_unwind(|| puruspe::invgammp(1e-300, 100.0)); // should be 0.03800698891694189
    println!("{r:?}");
}
```

Output (release and debug mode identical): `101.0`,
`2.0285042311847493e-28`, `1.881049407927448e-9`, `34.55874673932916`,
`2.5948963317332045e-20`, `0.6931471805547132`, `0.26505254963250424`,
`99.6668649180701`, then twice the panic message "Bad args in gammp" and
`Err(Any { .. })`.

## Expected

| (p, a) | Crate | mpmath (true value) | SciPy | Rel. error | Group |
|---|---|---|---|---|---|
| (1, 1) | 101.0 | inf (pole) | inf | inf | 1 |
| (0.1, 1000) | panic | 959.6939327288333 | 959.6939327288333 | - | 2 |
| (0.5, 0.01) | 2.0285042311847493e-28 | 4.465535018910355e-31 | 4.465535018910562e-31 | 453 | 3 |
| (1e-50, 2) | 1.881049407927448e-9 | 1.4142135623730952e-25 | 1.4142135623731e-25 | 1.3e16 | 4 |
| (1e-300, 100) | panic | 0.03800698891694189 | 0.038006988916941875 | - | 4 |
| (1 - 1e-15, 1) | 34.55874673932916 | 34.53957599234088 | 34.53957599234088 | 5.6e-4 | 5 |
| (1e-20, 1) | 2.5948963317332045e-20 | 1e-20 | 1.0000000000000022e-20 | 1.6 | 6 |
| (0.5, 1) | 0.6931471805547132 | 0.6931471805599453 | 0.6931471805599455 | 7.5e-12 | 6 |
| (0.99, 0.01) | 0.26505254963250424 | 0.26505255025158186 | 0.2650525502515818 | 2.3e-9 | 6 |
| (0.5, 100) | 99.6668649180701 | 99.66686491931549 | 99.66686491931549 | 1.2e-11 | 6 |

Second opinion: SciPy (`scipy.stats.gamma.ppf`) agrees with mpmath within
the tolerance at 391 of the 397 points; at the other 6 (p = 5e-324 for
a >= 100, and a = 1e6, p = 1e-10) SciPy is wrong too (verdict "both
wrong"), and the dataset's route is tested there by P(a, x) = p to 40
digits (see step 2). Closed forms: x = -ln(1 - p) for a = 1, and
x = (erfc^-1(1 - p))^2 for a = 1/2.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `puruspe = "=0.4.4"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The quantile is found by bisection and
   Newton's method on log P(a, e^u) = log p (log Q(a, e^u) = log(1 - p) for
   p > 1/2), on the cdf route of finding 012; tests check P(a, ppf(p)) = p
   (or Q = 1 - p) at a = 0.01 to 1e6 and p = 5e-324 to 1 - 1e-15, and the
   closed forms for a = 1 and 1/2 (`reference/tests/test_functions.py`,
   finding 013).
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points with finding 013's design and a new
   seed (Python `random.Random(20261104)`): a = 10^u with u uniform in
   [-2, 6], 6 significant digits; p uniform in (0, 1) (40%), p = 10^v with
   v uniform in [-300, -1] (30%), or p = 1 - 10^v with v uniform in
   [-15.5, -1] (30%); reference by the dataset's route at 30 digits.
   `invgammp` fails 2,717 (91%); every failure falls in a group above:
   - group 2: 640 panics after `gammp` returned NaN (a from 107 to
     9.8e5, p up to 0.494);
   - group 3: 101 (a = 0.0100 to 0.955, p up to 0.657: for a near 0.01
     even the median is below 1e-30), up to 452x;
   - group 4: 251 wrong answers (a = 1.02 to 561, p below 1.9e-21), up to
     5.7e243 relative, and 26 panics where the density at the start
     underflowed (a = 75 to 273, p below 1.3e-162);
   - group 5: 656 (p above 1 - 1.3e-5, every range of a), up to 9.3%;
   - group 6: 1,043 that pass when P is exact (below), up to 3.2e-9.
5. **Root cause** (from the source, `src/gamma.rs`, `invgammp`, a port of
   Numerical Recipes' `invgammp`), confirmed with copies (scratch, not
   committed):
   - A Rust copy of `invgammp` calling `puruspe::gammp` matches
     `puruspe::invgammp` bit for bit at all 432 grid points (the same
     values and the same panics). Instrumented, it shows the NaN first
     coming from `gammp` at the 61 group 2 points (and 640 of the 666
     random panics), and from a density of exactly 0 at the other 2 + 26.
   - A Python copy of the iteration (the same arithmetic in doubles) with
     P and Q from mpmath, correctly rounded, isolates invgammp's own code:
     it fails exactly the 107 points of groups 1, 3, 4 and 5 (plus one
     point, (1 - 1e-5, 0.5), where puruspe's inexact P happens to land
     within tolerance), so the 290 points of groups 2 and 6 are `gammp`'s.
     Switching on one remedy at a time fixes each group: +inf at p = 1
     (group 1); a relative stopping test |step| < x 2^-52 (group 3);
     more steps (200) and, where the 1e-3 clamp acts, the start
     x = (p Gamma(a + 1))^(1/a) from P ~ x^a / Gamma(a + 1) (group 4);
     solving Q(a, x) = 1 - p for p > 1/2 (group 5). Off the grid each of
     the 1,008 own-cause failures is fixed by exactly the remedy of its
     group (101, 251 and 656).
   - With exact P and Q and all remedies, 6 grid points still fail, all at
     p = 5e-324 with a >= 50, and 9 of the 3,000 random points (a = 116 to
     273, p below 2.4e-207). There P is subnormal or the density underflows
     to 0, so any step on P itself fails; only a step on log P (as the
     reference and SciPy's Boost route use) would work. SciPy also fails 5
     of these 6 grid points.

## Notes

The upper-tail start (Numerical Recipes' Wilson-Hilferty guess with
1 - p) is fine; only the equation solved is wrong there. `invgammp`
takes (p, a), the reverse of `gammp(a, x)`'s order. Finding 013 (statrs
`Gamma::inverse_cdf`) has an independent implementation with an absolute
stopping test as well. Boost.Math's `gamma_p_inv` starts from DiDonato
and Morris's approximation and uses Halley steps with a relative
tolerance, choosing P or Q by p.
