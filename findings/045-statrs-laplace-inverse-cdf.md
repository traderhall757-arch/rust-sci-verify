# 045: statrs Laplace inverse_cdf panics at p = 0 and 1, and returns inf where scale times the log overflows

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Laplace`,
  method `inverse_cdf` of `ContinuousCDF`)
- **Severity:** high (silently wrong answers: with location -2^1023 and
  scale 2^1023, `inverse_cdf(0.95)` is inf instead of 1.17e308), plus a
  panic at p = 0 and p = 1 for every location and scale
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 13 of 174 points of `reference/data/laplace_ppf.jsonl`
  (tolerance 1e-13). The exact list is `LAPLACE_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two causes:
  1. **Panic at p = 0 and p = 1** (12 grid points, tagged pole: every
     location and scale): `inverse_cdf` panics with the message
     "p must be in [0, 1]" although 0 and 1 are in that range and the
     quantiles there are -inf and +inf (statrs's own `Gumbel`, `Normal`
     and `Exp` return the infinity). Conventions: a panic is always crate
     wrong.
  2. **b L above the largest double** (1 grid point: location -2^1023,
     scale 2^1023, p = 0.95): `self.scale * (2. - 2. * p).ln()` overflows
     to -inf although the quantile mu - b L is finite, and the result is
     inf.
  Every other point passes (worst 2.9e-16); the 16 points where the true
  quantile is beyond the largest double are tagged overflow and statrs
  returns the right infinity there.

`Laplace::cdf` and `sf` fail for another reason: [finding 044](044-statrs-laplace-cdf.md).
Cause 2 is the same as cause 2 of [finding 043](043-statrs-gumbel-inverse-cdf.md) (`Gumbel`).

## Reproduction

The program of [finding 044](044-statrs-laplace-cdf.md), run in its own
fresh project, prints these rows too (lines 5 to 8 of its output: inf,
-inf, then a panic with the message "p must be in [0, 1]" for each of
p = 0 and p = 1, caught by `catch_unwind` as `Err(Any { .. })`).

## Expected

| mu, b, p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| 0, 1, 0 (and every other mu, b) | panic | -inf | -inf | panic |
| 0, 1, 1 (and every other mu, b) | panic | inf | inf | panic |
| -2^1023, 2^1023, 0.95 | inf | 1.1708241396246928e308 | inf | inf |
| 4.67591e307, 1.50284e308, 0.14745759143266557, off the grid | -inf | -1.367478051182768e308 | -inf | inf |

Second opinion: at p = 0 and 1 SciPy (`scipy.stats.laplace.ppf`) gives the
infinities. At the other rows SciPy multiplies by the scale the same way
and is wrong, so those rows were also computed by a second mpmath route:
mu + b log(2p) and mu - b log(2 - 2p) at 1200 bits, which agree with the
table to every printed digit. In `reference/tests/test_functions.py` the
generator's route is tested against the root of cdf(x) = p found by the
secant method at 1200 bits from the cdf (exp only, no logarithm), at
location -2^1023, scale 2^1023 among others.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program of finding 044, release and debug mode).
2. **Arbiter checked**: mpmath at 50 and 100 digits round to the same double
   (the generator enforces this); the route is tested as above. p = 0 and 1
   are tagged pole by the rule in `_laplace_ppf_tags`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random (mu, b, p) (Python `random.Random(20261022)`,
   after the 20,000 points of finding 044), with mu and b drawn as there;
   p = 10^w with w uniform in [-320, 0] for a quarter, uniform in (0, 1)
   for a quarter, 1 - 10^w with w uniform in [-16, 0] for a quarter and
   1/2 +- 10^w with w uniform in [-17, -1] for a quarter (a p that rounds
   to 0 or 1 is replaced by 1/2); mpmath at 50 and 100 digits (they agreed
   at all points). `inverse_cdf` fails 1, returning -inf where the quantile
   is -1.37e308 (cause 2; the last row above). A second scan
   (`random.Random(20261023)`, interleaved with finding 044's targeted
   points): 2,000 points with mu = -10^u and b = 10^v, u and v uniform in
   [307.5, 308.25], 6 significant digits, p uniform in (0, 1). It fails
   55, every one returning inf where the true quantile is finite (cause
   2). Cause 1 is at two fixed inputs and needs no scan: the source
   rejects p = 0 and 1 for every parameter.
5. **Root cause** (from the source, `src/distribution/laplace.rs`):
   `inverse_cdf` starts with `if p <= 0. || 1. <= p { panic!("p must be in
   [0, 1]") }` (cause 1), then computes `self.location + self.scale *
   (2. * p).ln()` or `self.location - self.scale * (2. - 2. * p).ln()`
   (cause 2 is the product overflowing). Confirmed with a copy that returns
   -inf and +inf at p = 0 and 1 and b (mu / b + L) where b L overflows:
   it passes every grid point and all 22,000 random points (worst 1.5e-15
   on the first scan, 1.6e-14 on the second).
