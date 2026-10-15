# 017: special Gamma::digamma aborts the process for large negative x, loses all digits next to its zero

- **Crate / version:** special 0.14.2 (`special::Gamma::digamma`, the
  digamma function psi(x) = Gamma'(x) / Gamma(x), called as `x.digamma()`;
  default features)
- **Severity:** high (the whole process aborts on a stack overflow, which
  `catch_unwind` cannot catch; next to the zero a silently wrong answer, with
  the wrong sign at the double nearest it)
- **Status:** approved (Claude, checklist complete, 2026-10-15)
- **Affected inputs:** 7 of 949 points of `reference/data/digamma.jsonl`
  (tolerance 4e-13). The exact list is `DIGAMMA` in
  `crates/check-special/src/expected.rs`; the report is
  [reports/special-0.14.2.md](../reports/special-0.14.2.md). Two groups:
  1. **Stack overflow** (2 points: x = -1000000.5 and -999999999999999.5):
     the call recurses once per unit from x up to 8, so it aborts the
     process once x is far enough below 0 for the stack. Where depends on the
     stack size and the build: on an 8 MiB main-thread stack x = -100000.5
     still returns and -200000.5 aborts in a debug build, -200000.5 returns
     and -500000.5 aborts in a release build (threads spawned with
     `std::thread` get 2 MiB by default, so less). It always aborts for
     x <= -2^53 (where x + 1 == x, so the recursion never ends) and for
     x = -inf.
  2. **Next to the zero x0 = 1.46163214496836234...** (5 points): relative
     error 7.3e-12 and 3.2e-11 at x0 -+ 1e-5, 1.8e-6 and 2.4e-6 at
     x0 -+ 1e-10, and 2.2 at the double nearest x0, where the crate returns
     +1.1e-16 for the true -9.2e-17. Off the grid (below) points fail up to
     2.1e-3 from x0.

Everything else passes, including every negative non-integer on the grid
down to -10000.5 and off it down to -2e4 (worst 9.1e-14).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add special@=0.14.2
use special::Gamma;

fn main() {
    // Next to the zero of digamma at 1.4616321449683623...
    for x in [
        1.4616321449683622_f64, // the double nearest the zero
        1.4616321448683622,     // zero - 1e-10
        1.4616321450683623,     // zero + 1e-10
        1.4616221449683622,     // zero - 1e-5
        1.4616421449683623,     // zero + 1e-5
        -100000.5,
    ] {
        println!("digamma({x:?}) = {:?}", x.digamma());
    }
    // Recurses once per unit below 8: overflows the stack and aborts.
    println!("digamma(-1000000.5) = {:?}", (-1000000.5_f64).digamma());
}
```

This is the program that was run, in release and in debug mode (same
output). It prints the "Crate" column below, then

```
thread 'main' (27881) has overflowed its stack
fatal runtime error: stack overflow, aborting
```

and exits with SIGABRT (status 134).

## Expected

| x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 1.4616321449683622 (nearest x0) | 1.1102230246251565e-16 | -9.241265521729427e-17 | -9.241265521729427e-17 | 2.2 (wrong sign) |
| 1.4616321448683622 (x0 - 1e-10) | -9.676714984863111e-11 | -9.676732496840103e-11 | -9.676732496840103e-11 | 1.8e-6 |
| 1.4616321450683623 (x0 + 1e-10) | 9.676737189323603e-11 | 9.676714013423533e-11 | 9.676714013423533e-11 | 2.4e-6 |
| 1.4616221449683622 (x0 - 1e-5) | -9.676766731137043e-6 | -9.67676673120742e-6 | -9.67676673120742e-6 | 7.3e-12 |
| 1.4616421449683623 (x0 + 1e-5) | 9.676678178083264e-6 | 9.676678178388793e-6 | 9.676678178388794e-6 | 3.2e-11 |
| -100000.5 | 11.512935464924329 | 11.512935464924395 | 11.512935464924396 | 5.7e-15 (passes) |
| -1000000.5 | process aborts | 13.815511557963816 | 13.815511557963816 | crash |

Second opinion: SciPy (`scipy.special.digamma`) agrees with mpmath within
1.8e-16 at every row. The zero itself: `mpmath.findroot(mpmath.digamma, 1.46)`
rounds to 1.4616321449683623 (the same double as 1.4616321449683622), tested
in `reference/tests/test_functions.py`.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `special = "=0.14.2"`
   (default features, as in the check crate: `no_std`, which takes `ln` from
   the `libm` crate); the program above. A second small program confirmed the
   abort at x = -inf and x = -1e300 (release build), and that NaN and +inf
   return NaN and inf.
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. `mpmath.digamma` is tested in
   `reference/tests/test_functions.py` against closed forms at integers and
   half-integers (DLMF 5.4.14, 5.4.15), the recurrence DLMF 5.5.2 from
   x + n for negative x and next to the zero (evaluated at 60 digits, where
   the cancellation leaves 40), the reflection formula DLMF 5.5.4 at large
   negative x, and the asymptotic series DLMF 5.11.2 at large x.
3. **Second opinion**: see "Expected".
4. **Off-grid**: a scratch program drew 12,000 points with xorshift64*
   (seed 20261015), in this order: 2,000 with x = x0 + u, u uniform in
   [-1e-2, 1e-2]; 2,000 with x = x0 +- 10^v, v uniform in [-16, -1]; 3,000
   with x = 10^v, v uniform in [-5, 8]; 3,000 with x = -10^v, v uniform in
   [-3, 4]; 2,000 with x = -10^v, v uniform in [4, 15] (digamma called only
   where x > -2e4, so as not to abort the scan; 39 points). Reference mpmath
   at 50 digits, tolerance 4e-13:
   - next to the zero: 134 of 2,000 (uniform) and 1,658 of 2,000 (log
     spaced) fail, worst 2.2 at the double nearest x0; the farthest failing
     point is 2.1e-3 from x0. SciPy is within tolerance at all of them.
   - positive x: 0 of 3,000 fail (worst 1.6e-13, at 1.4597, near x0).
   - negative x: 0 of 3,039 fail (worst 9.1e-14).
   The abort was located by calling the release and debug builds at
   x = -n - 0.5 for n = 1e4, 2e4, 3e4, 5e4, 1e5, 2e5, 5e5 on the main
   thread (results in "Affected inputs").
5. **Root cause** (from the source, `src/gamma.rs`):
   - **Stack overflow** (group 1). For p <= 8 the code is
     `return (p + 1.0).digamma() - p.recip();`, the recurrence
     psi(x) = psi(x + 1) - 1/x (DLMF 5.5.2) applied by recursion, one call
     per unit, with no reflection formula for negative x. The call is not
     in tail position (the subtraction follows), so each step keeps a stack
     frame, about 8 - x frames in all. For x <= -2^53, x + 1 rounds back to
     x and the recursion cannot end at any stack size; the same for -inf.
     Using the reflection formula psi(x) = psi(1 - x) - pi cot(pi x) for
     negative x (as Boost and SciPy do), with cot formed from the exact
     input, would bound the work.
   - **Cancellation next to the zero** (group 2). For x <= 8 the result is
     psi(x + 7) - sum 1/(x + k): about 2.13 minus about 2.13 near x0. A
     float64 emulation of the same steps (NumPy, scratch, not committed)
     reproduces the crate bit for bit at 3,976 of the 4,000 scan points
     next to the zero (the rest differ in the last bit of `ln`), and the
     crate's largest absolute error there is 1.04e-15, 2.3 ulps of 2.13:
     ordinary rounding, which is enormous relative to a result near 0. Boost
     (`digamma_imp_1_2`) and SciPy expand around x0 in x - x0 with x0 held
     to more than double precision, which keeps the relative error to a few
     ulps; special has no such branch.

The check runs `digamma` with `Check::run_isolated` (child processes), so
the abort is recorded as "crash: stack overflow" instead of ending the
test.
