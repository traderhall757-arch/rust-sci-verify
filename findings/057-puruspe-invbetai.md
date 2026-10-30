# 057: puruspe invbetai panics where the density underflows, starts far from tail quantiles, solves the upper tail through P, and inherits betai's errors

- **Crate / version:** puruspe 0.4.4 (`puruspe::invbetai(p, a, b)`, the x
  with I_x(a, b) = p: the inverse of the regularized incomplete beta
  function, i.e. the quantile of the beta distribution)
- **Severity:** high (silently wrong answers: `invbetai(1e-300, 1, 1)` is
  0 instead of 1e-300, `invbetai(1 - 1e-15, 1, 1000)` is 0.393 instead of
  0.0339; a panic "Bad x in routine betai" for ordinary inputs such as
  `invbetai(1e-20, 1000, 1)` and `invbetai(0.5, 1e4, 1e4)`)
- **Status:** approved (Claude, checklist complete, 2026-10-30)
- **Affected inputs:** 1035 of the 1620 points of
  `reference/data/betaincinv.jsonl` (tolerance 2e-13); every one of the 81
  pairs (a, b) on the grid fails somewhere. The exact lists are
  `INVBETAI` (244 points, groups 1 to 4: invbetai's own code) and
  `INVBETAI_BETAI` (791 points, groups 5 and 6: where an error of the
  `betai` it calls shows first; 3 of them would fail even with an exact P,
  see step 5) in `crates/check-puruspe/src/expected.rs`;
  the report is [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md).
  Six groups, by cause (checklist step 5):
  1. **Panic where the density underflows at the start** (42 points,
     a = 5 to 1e4, b = 1 to 100, p <= 1e-20; 17 of them at p = 5e-324):
     for a, b >= 1 the start is Numerical Recipes' normal-approximation
     guess, which in the far lower tail is far below the root (6.5e-102
     for (5e-324, 5, 1), true 2.2e-65). The density
     exp((a - 1) ln x + (b - 1) ln(1 - x) + afac) then underflows to
     exactly 0, the step err / 0 is -inf, the damped Halley step
     -inf / +inf is NaN, and the next `betai(a, b, NaN)` fails its
     argument check and panics with "Bad x in routine betai"
     (`invbetai(1e-20, 1000, 1)`, true p^(1/a) = 0.955).
  2. **A start far from a tail quantile and 10 steps** (66 points,
     a, b >= 1): the same guess is far from the root in the far lower tail
     (59 points, p <= 1e-10) and, for a = 1 and b >= 5, in the far upper
     tail (7 points, p >= 1 - 1e-10). 57 use all 10 steps without meeting
     the stopping test: `invbetai(1e-100, 10, 1)` = 4.9e-22, true
     p^(1/a) = 1e-10; `invbetai(1 - 1e-15, 1, 1000)` = 0.393, true
     1 - (1 - p)^(1/b) = 0.0339. At 7 points (a = 1, p = 1e-300) the guess
     itself rounds to 0 (e^(2w) overflows) and the loop returns it
     unchanged: `invbetai(1e-300, 1, 1)` = 0, true 1e-300 (the uniform
     distribution). The other 2, (5e-324, 2, 10) and (1 - 1e-15, 1, 5),
     meet the stopping test while still wrong. Relative errors 5.8e-7 to
     2.3e8. The 7 upper-tail points also need group 3's remedy: for a = 1
     the start remedy (below) is the exact answer, so it fixes them without
     a step, but 200 steps from the old start fail all 7, and group 3's
     remedy with 200 steps fixes 6.
  3. **Upper tail through P** (105 points, b >= 2, p >= 1 - 1e-5): the
     equation solved is I_x(a, b) - p = 0. Near p = 1 the doubles are
     2^-53 apart, so even a correctly rounded P leaves x uncertain by
     about 2^-54 / f(x), f the density; it should solve
     1 - I_x(a, b) = 1 - p (exact for p >= 1/2) instead. Errors 2.2e-13 to
     17 (`invbetai(1 - 1e-15, 0.5, 100)` = 0.27585, true 0.27599). 5
     points (a = 1 and 2 with b >= 100) also need group 2's remedy. At
     p = 1 - 1e-5 this floor is close to the tolerance, so there the split
     between groups 2, 3 and 6 is a matter of which error is larger: with
     a correctly rounded log-gamma in `afac` instead of `math.lgamma`, the
     exact-P copy moves (1 - 1e-5, 2, 1e4) to passing and
     (1 - 1e-5, 1, 1e4) to failing (groups 2 and 3 would be 68 and 102).
     The lists are not affected: the harness checks their union.
  4. **p = 5e-324** (31 points, a >= 2): p is the smallest subnormal, and
     the representable values of P near it are 4.9e-324 apart, so the
     step's err = P - p is known only to about 100%, and x only to about
     1/a: `invbetai(5e-324, 1000, 1000)` = 0.138493, true 0.138438
     (3.9e-4). Errors 9.8e-6 (a = 1e4) to 1.5e9 (a = 5, where the start
     is also far from the root, group 2's cause).
  5. **Where finding 001's betai bugs show** ([finding 001](001-puruspe-betai.md),
     93 points, 12 of them panics): the iteration solves for puruspe's P.
     A point is in this group if a = b = 1e4 or finding 001's fixed
     `betai` (`crates/puruspe-fix`) makes it pass or cuts its error more
     than tenfold. At a = b = 1e4 (18 points) `betai` returns 0 or -0 for
     every x (bug 1): 9 points panic as in group 1 after the resulting
     huge step (`invbetai(0.5, 1e4, 1e4)`, true 0.5), the rest are up to
     4.5% off; with the fixed `betai`, (5e-324, 1e4, 1e4) keeps its 5.7e-4
     error (group 4's cause) and (1 - 1e-15, 1e4, 1e4) still misses by
     5.1e-6 (group 3's cause). At 75 points with a or b <= 3000 the
     continued-fraction switch (bug 3, wider scope) shows: 24 crash or are
     off by more than 1e-8 (`invbetai(1e-20, 1000, 1000)` = 0, true 0.398;
     3 panics at a = 1e4 with b = 100 and 1000), while 22 are off by less
     than 1e-11, the size of group 6, and 16 pass when only `betai`'s
     prefactor is made exact (so for those, bugs 3 and 4 overlap). With the
     fixed `betai` 52 of the 93 pass.
  6. **betai's accuracy** (698 points, 62 of the 81 pairs): the iteration
     solves for puruspe's P, and with an exact P all 698 pass. 696 lose
     2.1e-13 to 2.3e-10: `betai(2, 2, 0.5)` = 0.4999999999976814, so
     `invbetai(0.5, 2, 2)` = 0.5000000000015458 (true 0.5);
     `betai(0.1, 1, 1e-200)` = 1e-20 (1 + 1.2e-11), so
     `invbetai(p, 0.1, 1)` = p^10 (1 - 1.2e-10) for every p <= 0.9. With
     `betai`'s prefactor made exact (finding 001, bug 4: built from
     puruspe's `ln_gamma`, [finding 033](033-puruspe-ln-gamma.md)) 663 of
     the 696 pass; the other 33, 27 of them at b = 1e4, still miss by
     2e-13 to 6.9e-12 (e.g. (0.1, 0.1, 1e4) at 3.3e-12, as (p, a, b)) from other
     rounding inside `betai`, not traced further (likely b ln(1 - x)
     with 1 - x rounded). (The exact prefactor here is the double nearest
     the mpmath value; `betai` with a, b > 3000 is left as is.) At 420 of the 698 the loop also uses all 10
     steps without meeting its stopping test. The last 2 points (a = 1,
     b = 10 and 100, p = 1 - 1e-5; 2.8e-13 and 3.9e-13) are at group 3's
     floor: every x within 2.6e-13 (4.6e-13 for b = 100) of the root gives
     the same rounded P, wider than the tolerance, so where the loop
     settles depends on its path (puruspe's P equals p exactly at the
     returned x).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // invbetai(p, a, b): the x with I_x(a, b) = p
    let r = std::panic::catch_unwind(|| puruspe::invbetai(1e-20, 1000.0, 1.0)); // should be 0.954992586021436 (= p^(1/a))
    println!("{r:?}");
    println!("{:?}", puruspe::invbetai(1e-300, 1.0, 1.0)); // should be 1e-300 (= p)
    println!("{:?}", puruspe::invbetai(1e-100, 10.0, 1.0)); // should be 1e-10 (= p^(1/a))
    println!("{:?}", puruspe::invbetai(0.999999999999999, 1.0, 1000.0)); // should be 0.03394989346167808 (= 1 - (1 - p)^(1/b))
    println!("{:?}", puruspe::invbetai(0.999999999999999, 0.5, 100.0)); // should be 0.275986209102702
    println!("{:?}", puruspe::invbetai(5e-324, 1000.0, 1000.0)); // should be 0.1384383837250825
    let r = std::panic::catch_unwind(|| puruspe::invbetai(0.5, 10000.0, 10000.0)); // should be 0.5
    println!("{r:?}");
    println!("{:?}", puruspe::invbetai(1e-20, 1000.0, 1000.0)); // should be 0.3975320608971773
    println!("{:?}", puruspe::invbetai(0.5, 2.0, 2.0)); // should be 0.5
}
```

Output (release and debug mode identical): the panic message "Bad x in
routine betai" and `Err(Any { .. })`, `0.0`, `4.868566054987597e-22`,
`0.39342642972360853`, `0.2758533347864622`, `0.1384928315311716`, again
the panic message and `Err(Any { .. })`, `0.0`, `0.5000000000015458`.

## Expected

| (p, a, b) | Crate | mpmath (true value) | SciPy | Rel. error | Group |
|---|---|---|---|---|---|
| (1e-20, 1000, 1) | panic | 0.954992586021436 | 0.954992586021436 | - | 1 |
| (1e-300, 1, 1) | 0.0 | 1e-300 | 1e-300 | 1 | 2 |
| (1e-100, 10, 1) | 4.868566054987597e-22 | 1e-10 | 9.999999999999987e-11 | 1.0 | 2 |
| (1 - 1e-15, 1, 1000) | 0.39342642972360853 | 0.03394989346167808 | 0.03394989346167808 | 10.6 | 2 |
| (1 - 1e-15, 0.5, 100) | 0.2758533347864622 | 0.275986209102702 | 0.27598620910270205 | 4.8e-4 | 3 |
| (5e-324, 1000, 1000) | 0.1384928315311716 | 0.1384383837250825 | 0.1384383922415156 | 3.9e-4 | 4 |
| (0.5, 1e4, 1e4) | panic | 0.5 | 0.5 | - | 5 |
| (1e-20, 1000, 1000) | 0.0 | 0.3975320608971773 | 0.39753206089717746 | 1 | 5 |
| (0.5, 2, 2) | 0.5000000000015458 | 0.5 | 0.5 | 3.1e-12 | 6 |

Second opinion: SciPy (`scipy.special.betaincinv`) agrees with mpmath
within the tolerance at 971 of the 1035 points; at the other 64 SciPy is
wrong too (verdict "both wrong"): 31 of group 4 (p = 5e-324, where SciPy
also fails), 15 of group 1, 8 of group 2, 9 of group 6 and 1 of group 5,
62 of them with p <= 1e-300 or a quantile below the smallest normal
double, and SciPy's 2.8e-13 and 2.1e-13 at (a, b) = (1000, 1e4),
p = 0.999 and 1 - 1e-5 (the dataset documents SciPy's 105 failures). Closed forms: x = p^(1/a) for
b = 1, x = 1 - (1 - p)^(1/b) for a = 1, x = 1/2 at p = 1/2 for a = b.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `puruspe = "=0.4.4"`
   (the program above, release and debug mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The quantile is found by Newton's method on
   log I in s = log(x / (1 - x)) (log(1 - I) for p > 1/2), with the
   continued fraction on whichever side converges; tests check it against
   the series and the closed forms for b = 1, a = 1 and a = b = 1/2
   (`reference/tests/test_functions.py`, finding 015).
3. **Second opinion**: see "Expected". Where SciPy is wrong too, the
   table's row (5e-324, 1000, 1000) was checked by a second route,
   `mpmath.betainc` at 40 digits: I at the reference x is
   p (1 + 6.8e-14), and the next double up gives p (1 + 2.4e-13), so the
   reference is the closest double.
4. **Off-grid**: 3,000 random points (Python `random.Random(20261105)`):
   a and b = 10^u with u uniform in [-1, 4], 6 significant digits; p
   uniform in (0, 1) (40%), p = 10^v with v uniform in [-300, -1] (30%),
   or p = 1 - 10^v with v uniform in [-15.5, -1] (30%); reference by the
   dataset's route at 30 digits. `invbetai` fails 2,158 (72%); every
   failure falls in a group above:
   - group 1: 89 panics (a = 11.7 to 9,990, b = 1.04 to 33, p from
     1.2e-300 to 2.2e-48);
   - group 2: 246 (a = 0.56 to 9,370, b = 0.70 to 9,440, p below
     2.3e-21), up to 5.9e7 relative; 2 of them fail even with every
     remedy (step 5);
   - group 3: 437 (p above 1 - 1.4e-5), up to 4.2; 424 fixed by group 3's
     remedy alone, 13 also need group 2's;
   - group 4: none (the scan's smallest p is about 1e-300);
   - group 5: 153, 22 of them panics;
   - group 6: 1,233 that pass when P is exact, up to 1.4e-10 (1 of them,
     p near 1, is at group 3's floor like group 6's last 2 grid points).
5. **Root cause** (from the source, `src/beta.rs`, `invbetai`, a port of
   Numerical Recipes' `invbetai`), confirmed with copies (scratch, not
   committed):
   - A Rust copy of `invbetai` calling `puruspe::betai` matches
     `puruspe::invbetai` bit for bit at all 1620 grid points and all
     3,000 random points (the same values and the same panics; debug and
     release builds identical). Instrumented, it shows the density
     exactly 0 at the start at all 42 group 1 points, and the step then
     NaN. Swapping in finding 001's fixed `betai` separates group 5.
     Swapping `libm::lgamma` into `afac` makes group 6's last 2 points
     pass, but only by changing the path (an error in `afac` scales the
     steps, not the root); see group 6.
   - A Python copy of the iteration in doubles, with P and 1 - P from
     mpmath (correctly rounded) and Python's `math.lgamma` (accurate to a
     few units in the last place, not correctly rounded) in `afac`,
     isolates invbetai's own code. It fails 248 grid points: the 244 of
     groups 1 to 4, 3 of group 5 where betai's error comes first
     (a = b = 1e4 at p = 5e-324, 1 - 1e-10 and 1 - 1e-15), and 1,
     (1 - 1e-5, 0.5, 100), where puruspe's inexact P happens to land
     within tolerance. Every other
     point of groups 5 and 6 passes, so those failures are `betai`'s.
   - Remedies switched on one at a time fix the groups: solving
     1 - I = 1 - p for p > 1/2 (group 3; 5 also need the next one); a
     start at or above the lower bound (p a B(a, b))^(1/a), from
     I_x(a, b) <= x^a / (a B(a, b)) for b >= 1, mirrored in the upper tail
     for a >= 1, with 200 steps (groups 1 and 2); at p = 5e-324 also a
     step on log P computed without rounding P to a double (group 4, and
     11 of group 1's 17 points at p = 5e-324). With all of them the copy
     passes all 1620 grid points. Off the grid all but 2 of the 3,000
     pass: two quantiles below the smallest normal double (8.4e-315 and
     6.8e-311, a < 1), which still miss.

## Notes

For a < 1 or b < 1 the start is Numerical Recipes' power-law guess, and no
grid failure comes from it (groups 1 and 2 have a, b >= 1 on the grid).
Off the grid 3 of group 2's points use it; all 3 have quantiles below the
smallest normal double, and 2 of them are the two that no remedy fixes.
The stopping test |step| < x 2^-52 is relative and fine. statrs's
`Beta::inverse_cdf` on the same data is finding 020; its causes differ
(it loops forever; puruspe stops after 10 steps).
