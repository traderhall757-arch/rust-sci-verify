# 018: special Gamma::trigamma is finite at its poles and inaccurate for negative x

- **Crate / version:** special 0.14.2 (`special::Gamma::trigamma`, the
  trigamma function psi'(x), called as `x.trigamma()` (it takes `&self`);
  default features)
- **Severity:** high (wrong answer silently: finite values from 8.5e28 to
  6.4e36 at the poles, and up to 90 times the true value at large negative x)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 31 of 943 points of `reference/data/trigamma.jsonl`
  (tolerance 1e-13). The exact list is `TRIGAMMA` in
  `crates/check-special/src/expected.rs`; the report is
  [reports/special-0.14.2.md](../reports/special-0.14.2.md). Two groups,
  one cause:
  1. **Finite at the poles** (30 points): every negative integer from -1 to
     -30 gives a finite number (6.6e32 at -1, 8.5e28 at -30) where trigamma
     has a pole; a finite value at a pole is crate wrong
     (`docs/conventions.md`). x = 0 gives inf (accepted).
  2. **Accuracy loss for negative non-integers** (1 grid point, 4.8e-13 at
     x = -1000.25): the error grows with |x|. Off the grid (below) it starts
     near x = -1 and reaches 8.1e-10 at -6895, 1.4e-7 below -1e6, 0.2 (20%)
     below -1e11, and 90 times the true value near -9.5e14.

Positive x passes everywhere (worst 2.1e-15 on and off the grid).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add special@=0.14.2
use special::Gamma;

fn main() {
    // Poles: the true value is infinite.
    for x in [-1.0_f64, -2.0, -30.0] {
        println!("trigamma({x:?}) = {:?}", x.trigamma());
    }
    for x in [-1000.25_f64, -6894.998082352342, -1000000.25, -948584831202026.1] {
        println!("trigamma({x:?}) = {:?}", x.trigamma());
    }
}
```

This is the program that was run, in release and in debug mode (same
output); it prints the "Crate" column below.

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| -1 | 6.580790147320947e32 | pole (+inf) | inf | finite at a pole |
| -2 | 1.6451975368302367e32 | pole (+inf) | inf | finite at a pole |
| -30 | 8.494017819538807e28 | pole (+inf) | inf | finite at a pole |
| -1000.25 | 19.738209551709215 | 19.738209551699786 | 19.738209551699782 | 4.8e-13 |
| -6894.998082352342 | 271936.5758642824 | 271936.576084019 | 271936.57608401886 | 8.1e-10 |
| -1000000.25 | 19.739207812190273 | 19.739207802179468 | 19.739207802179465 | 5.1e-10 |
| -948584831202026.1 | 6130.72057765493 | 67.39387440291296 | (does not return) | 90 |

Second opinion: SciPy (`scipy.special.polygamma(1, x)`) agrees with mpmath
within 4.3e-16 where it returns. At x = -948584831202026.1 it does not
return (its Hurwitz zeta sums one term per unit from x up), so the second
route is the reflection formula itself: this double is
-948584831202026.125, so sin^2(pi x) = sin^2(pi / 8) and
pi^2 / sin^2(pi / 8) = 67.3938744 (the other term, psi'(1 - x), is about
1e-15).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `special = "=0.14.2"`
   (default features, as in the check crate: `no_std`, which takes `sin` and
   `pow` from the `libm` crate); the program above.
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The reference is `mpmath.psi(1, x)` for x > 0
   and the reflection formula with `mpmath.sinpi` for x < 0, because
   `mpmath.psi(1, x)` recurs one step at a time up from negative x (5 s at
   -1e6, endless at -1e15). Tested in `reference/tests/test_functions.py`
   against closed forms at integers and half-integers (DLMF 5.15.2, 5.15.5),
   the recurrence from x + n for negative x, `mpmath.psi(1, x)` itself for
   x from -0.25 to -10000.5, and the asymptotic series DLMF 5.15.8 at large
   x. Poles are tagged by rule (non-positive integers) and not evaluated.
3. **Second opinion**: see "Expected".
4. **Off-grid**: the scan of [finding 017](017-special-digamma.md) (12,000
   points, xorshift64*, seed 20261015; same ranges), with trigamma at every
   point; reference mpmath at 50 digits, tolerance 1e-13:
   - positive x: 0 of 7,000 fail (worst 2.1e-15).
   - x = -10^v, v in [-3, 4]: 589 of 3,000 fail, by decade of |x|: 0 below
     0.1, 1 of 417 in [0.1, 1) (the failing point closest to 0 is -0.998),
     3 of 414 in [1, 10) (worst 8.8e-13), 32 of 421 in [10, 100), 192 of
     429 in [100, 1e3) (7.7e-11), 361 of 437 in [1e3, 1e4) (8.1e-10).
   - x = -10^v, v in [4, 15]: 1,985 of 1,987 non-integers fail: worst 4.1e-9
     in [1e4, 1e5), 1.4e-7 in [1e6, 1e7), 1.3e-3 in [1e10, 1e11), 0.20 in
     [1e11, 1e12), 90 in [1e14, 1e15).
   SciPy is within tolerance at every failing point of the first range (it
   was not run below -1e4).
5. **Root cause** (confirmed with a copy of `trigamma`, scratch, not
   committed, that uses the same `libm` functions, returns the crate's
   result bit for bit at all 5,155 negative points of the grid and the
   scan, and can switch the cause off). For x <= 0 the code uses the
   reflection formula psi'(x) = (pi / sin(pi x))^2 - psi'(1 - x), with
   `(PI * x).sin()`. The product pi x is rounded before `sin` sees it: an
   absolute error up to |x| pi 2^-53 in the angle, which the steep sin^2
   near its zeros turns into the growing relative error of group 2. At a
   negative integer n, pi n rounds to a double that is not a multiple of pi,
   so sin returns about 1e-16 |n| instead of 0, and the result is finite
   (group 1). Forming sin(pi x) from the exact input (reduce r = x - 2
   round(x / 2), exact, then r to sign(r) - r if |r| > 1/2, then
   sin(pi r), as a `sinpi` does) removes every failure: 0 of 5,111 negative
   non-integers fail (worst 6.4e-16), and every pole gives +inf.
