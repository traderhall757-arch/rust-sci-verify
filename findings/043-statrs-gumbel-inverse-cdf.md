# 043: statrs Gumbel inverse_cdf loses all digits next to p = 1/e, and returns inf where scale times log(-log p) overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Gumbel`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location 0 and scale 1,
  `inverse_cdf(0.36787944117144233)` is 0 instead of 3.4e-17; with location
  -2^1023 and scale 2^1023, `inverse_cdf(0.9)` is inf instead of 1.12e308),
  plus an accuracy loss of about 1.1e-16 / |e p - 1| next to p = 1/e for
  ordinary parameters
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 20 of 192 non-pole points of
  `reference/data/gumbel_ppf.jsonl` (204 points, tolerance 1e-13). The exact
  list is `GUMBEL_PPF` in `crates/check-statrs/src/expected.rs`, the report
  is [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes:
  1. **p next to 1/e** (18 grid points: location 0, every scale, p = the
     doubles nearest 1/e and 1/e +- 10^-k for k = 5, 10, 15, except
     1/e - 1e-5, which passes by luck of rounding): the quantile is
     mu - beta L with L = log(-log p), and -log p is close to 1 there, so
     the log subtracts two nearly equal numbers and the absolute error of
     about 1.1e-16 in -log p becomes a relative error of about
     1.1e-16 / |e p - 1| in L. Where L is the quantile's main part (mu = 0,
     or beta |L| not small next to |mu|) that is the result's error: 0
     instead of 3.4e-17 at the double nearest 1/e, and with location 0 and
     scale 10, `inverse_cdf(0.3679)` is 3.8e-13 off. The point p = 1/e is
     the mode's quantile (37% of the mass lies below it), so ordinary
     inputs are hit.
  2. **beta L above the largest double** (2 grid points: location -2^1023,
     scale 2^1023, p = 0.9 and 0.95): `self.scale * L` overflows to -inf
     although mu - beta L is finite, and the result is inf.
  Every other non-pole grid point passes (worst 4.5e-15); the 15 points
  where the true quantile is beyond the largest double are tagged overflow
  and statrs returns the right infinity there.

`Gumbel::cdf` and `sf` fail for other reasons: [finding 042](042-statrs-gumbel-cdf.md).

## Reproduction

The program of [finding 042](042-statrs-gumbel-cdf.md), run in its own
fresh project, prints these rows too (lines 3 to 6 of its output).

## Expected

| mu, beta, p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 0, 1, 0.36787944117144233 (nearest 1/e) | 0.0 | 3.3784855259134224e-17 | 0.0 | 1 |
| 0, 1, 0.36787944127144234 (1/e + 1e-10) | 2.718282266314957e-10 | 2.7182823912192446e-10 | 2.718282266314957e-10 | 4.6e-8 |
| 0, 10, 0.3679, off the grid | 0.000558846901118221 | 0.0005588469011184362 | 0.000558846901118221 | 3.8e-13 |
| -2^1023, 2^1023, 0.9 | inf | 1.1238883801828629e308 | inf | inf |
| -2^1023, 2^1023, 0.95 | inf | 1.7709032367707242e308 | inf | inf |

Second opinion: SciPy (`scipy.stats.gumbel_r.ppf`) computes
-log(-log p) the same way and is wrong at every row, so every row was also
computed by a second mpmath route: mu - beta log(-log p) at 1200 bits,
which agrees with the table to every printed digit. In
`reference/tests/test_functions.py` the generator's route is tested against
the root of cdf(x) = p found by the secant method at 1200 bits from the cdf
(exp only, no logarithm), at these parameters and at the 1/e probes; a
second test checks that the cancellation next to 1/e stays within the
generator's 60 extra bits.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program of finding 042, release and debug mode).
2. **Arbiter checked**: mpmath at 50 and 100 digits round to the same double
   (the generator enforces this); the route is tested as above.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, beta, p) (Python
   `random.Random(20261020)`, after the 20,000 points of finding 042), with
   mu and beta drawn as there; p = 10^w with w uniform in [-320, 0] for a
   quarter, uniform in (0, 1) for a quarter, 1 - 10^w with w uniform in
   [-16, 0] for a quarter and (1 +- 10^w) / e with w uniform in [-17, -2]
   for a quarter (a p that rounds to 0 or 1 is replaced by 1/2); mpmath at
   50 and 100 digits (they agreed at all points). `inverse_cdf` fails
   2,762: 2,758 with |e p - 1| <= 1.0e-3 (cause 1, worst error 1.0) and 4
   returning inf (cause 2). A second scan (`random.Random(20261021)`, after
   the targeted points of finding 042): 2,000 points with |mu| = 10^u and
   beta = 10^v, u and v uniform in [307.5, 308.25], p uniform in (0, 1).
   It fails 59: 57 returning inf where the true quantile is 1.2e307 to
   1.78e308 (cause 2), and 2 where mu - beta L itself cancels (mu and
   beta L agree to about 1 part in 1,700; see Notes).
5. **Root cause** (from the source, `src/distribution/gumbel.rs`):
   `inverse_cdf` computes `self.location - self.scale * ((-(p.ln())).ln())`.
   Cause 1 is the outer `ln` of a number close to 1; cause 2 is the product
   `self.scale * L` overflowing. Confirmed with a copy that, for
   0.75 <= e p <= 1.25, forms d = e p - 1 with e as a double-double (one
   fused multiply-add; e p - 1 is then exact up to e's low part) and
   L = ln_1p(-ln_1p(d)); and that returns beta (mu / beta - L) where
   beta L overflows. This copy passes every grid point (worst 4.9e-16) and
   all 20,000 random points (worst 9.6e-15).

## Notes

The copy also fails the 2 points of the second scan where mu - beta L
cancels: any method that rounds L to a double loses (|mu| + |beta L|) /
|mu - beta L| ulps there (about 1,700, so 1.9e-13), and statrs gives the
same value as the copy. Those 2 points are not counted as part of this
finding. On the grid that ratio is at most about 14 (tested), within the
tolerance; a future off-grid failure from it alone would not by itself be
a bug.
