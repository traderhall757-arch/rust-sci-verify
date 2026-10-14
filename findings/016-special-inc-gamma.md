# 016: special Gamma::inc_gamma returns 1 for every x above 1e8, loses digits for large shape

- **Crate / version:** special 0.14.2 (`special::Gamma::inc_gamma`, the
  regularized lower incomplete gamma function P(a, x), called as
  `x.inc_gamma(a)`; default features)
- **Severity:** high (wrong answer silently: exactly 1.0 for every x > 1e8,
  where the true value can be 0 or anything in between)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 46 of 781 points (tolerance 1e-12): 1 of 192 in
  `reference/data/gammainc.jsonl`, 19 of 538 in `gamma_cdf.jsonl`, 26 of 51
  in `gammainc_large.jsonl`. The exact lists are `INC_GAMMA`,
  `INC_GAMMA_CDF` and `INC_GAMMA_LARGE` in
  `crates/check-special/src/expected.rs`; the report is
  [reports/special-0.14.2.md](../reports/special-0.14.2.md). Three groups,
  by cause (section "Root cause"):
  1. **x > 1e8 returns 1** (21 points, all in `gammainc_large`): wrong
     whenever the shape is above about 1e8 - 8.5e4 (so that P(a, x) is not 1
     to double precision). Errors are total: 1 instead of 0 at
     x = a / 2, 1 instead of 0.5 at x = a.
  2. **Digits lost for large shape** (23 points, a >= 1e4 on the grid,
     a >= 800 off it): the error grows like the shape, worst 1.1e-11
     (a = 1e4), 1.9e-10 (1e5), 1.2e-9 (1e6) and 2.0e-7 (1e8, at
     x = a - sqrt(a), where x is not yet above 1e8).
  3. **Subnormal results flushed to 0** (2 points): P(200, 2) = 2.8e-316 and
     P(50, 1e-5) = 3.3e-315. Every result below 2.2e-308 is returned as 0.

`special::Gamma::gamma`, `ln_gamma` and the sign `ln_gamma` returns pass
every point of `gamma`, `gammaln` and `gammasgn` (they are libm's `tgamma`
and `lgamma_r`); see the report.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add special@=0.14.2
use special::Gamma;

fn main() {
    // P(a, x) is `x.inc_gamma(a)` in special.
    for (a, x) in [
        (1e9, 5e8),          // true P underflows to 0
        (1e9, 1e9),          // true P = 0.50000420522087
        (1e8, 1.0001e8),     // true P = 0.84134474647180
        (1e6, 997000.0),     // true P = 1.3381041673135997e-3
        (1e6, 1e6),          // true P = 0.5001329807608725
        (200.0, 2.0),        // true P = 2.7852507e-316 (subnormal)
        (50.0, 1e-5),        // true P = 3.28791718e-315 (subnormal)
    ] {
        println!("P({a:e}, {x:e}) = {:?}", x.inc_gamma(a));
    }
}
```

This is the program that was run, in release and in debug mode (same
output); it prints the "Crate" column below.

## Expected

| a | x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| 1e9 | 5e8 | 1.0 | 0 (about e^-1.5e8) | 0 | total |
| 1e9 | 1e9 | 1.0 | 0.5000042052208701 | 0.5000042052208701 | 1.0 |
| 1e8 | 1.0001e8 | 1.0 | 0.8413447464717988 | 0.8413447464717987 | 0.19 |
| 1e6 | 997000 | 0.0013381041688938313 | 0.0013381041673135997 | 0.0013381041673135986 | 1.2e-9 |
| 1e6 | 1e6 | 0.5001329811027386 | 0.5001329807608725 | 0.5001329807608725 | 6.8e-10 |
| 200 | 2 | 0 | 2.7852507e-316 | 0 | 1.3e-8 (relative to 2.2e-308) |
| 50 | 1e-5 | 0 | 3.28791718e-315 | 0 | 1.5e-7 (relative to 2.2e-308) |

Second opinion: SciPy (`scipy.special.gammainc`) agrees with mpmath within
1.1e-15 at the first five rows. It also flushes the two subnormal results
to 0; for those the finite sum Q(200, 2) = e^-2 sum_{k<200} 2^k / k! (DLMF
8.4.10, at 400 digits) gives P = 2.78525070506e-316, and
`mpmath.gammainc` and the first two series terms give 3.28791718e-315 for
P(50, 1e-5). Errors below 2.2e-308 are measured relative to that number
(`docs/conventions.md`).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `special = "=0.14.2"`
   (default features, as in the check crate: `no_std`, which takes `ln`,
   `exp` and `lgamma` from the `libm` crate); the program above.
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. `gammainc` uses `mpmath.gammainc`; `gamma_cdf`
   and `gammainc_large` use the series DLMF 8.7.1 and the continued
   fraction DLMF 8.9.2, tested in `reference/tests/test_functions.py`
   against `mpmath.gammainc`, the finite sum DLMF 8.4.10, closed forms, and
   (new) numerical integration of DLMF 8.2.1 / 8.2.2 at the bulk points of
   `gammainc_large`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: three random scans (Python `random.Random(16)`, in this
   order; the crate run from a scratch binary, mpmath at 30 digits):
   - **Shapes 1e8 to 1e12** (300 points): a = 10^u, u uniform in [8, 12],
     x = a + z sqrt(a), z uniform in [-8, 8]; reference by numerical
     integration (the tested route above). The crate returned exactly 1.0
     at all 300; the true P was below 0.5 at 147. 281 fail (the other 19
     have P within 1e-12 of 1), smallest failing a = 1.017e8, worst error
     5.7e14 (P = 1.7e-15 returned as 1).
   - **Shapes 10 to 1e6** (1,895 points; 2,000 drawn, the ones with x > 0
     kept): a = 10^u, u uniform in [1, 6], x = a + z sqrt(a), z uniform in
     [-12, 12]; reference `_mp_gamma_cdf`. 585 fail: 0 of 290 with
     a < 100 (worst 5.8e-14), 2 of 441 with 100 <= a < 1000 (worst 1.4e-12,
     smallest failing a = 800.8), 121 of 366 below 1e4 (worst 2.6e-11), 208
     of 377 below 1e5 (2.2e-10), 254 of 421 below 1e6 (worst 3.3e-9 at
     a = 781367, x = 781198).
   - **Subnormal results** (300 points): a uniform in [1, 500], x chosen so
     that x^a / Gamma(a + 1) is about 10^-u, u uniform in [310, 322]. The
     crate returned 0 at all 300; the reference was subnormal at 215; 141
     fail (the rest are below 2.2e-320, within tolerance of 0 under the
     metric), worst 2.6e-3.
5. **Root cause** (confirmed with a copy of `inc_gamma`, scratch, not
   committed, that uses the same `libm` functions, returns the crate's
   result bit for bit at all 781 points, and can switch each cause off).
   `inc_gamma` is AS 239 (Shea 1988), via John Burkardt's C version: a
   series for x <= 1 or x < p, a continued fraction otherwise, both
   multiplied by exp of a log prefix.
   - **The x > 1e8 shortcut** (group 1). The code returns 1 as soon as
     x > XBIG = 1e8. In AS 239 that test comes *after* the branch for
     p > PLIMIT = 1000 (a Wilson-Hilferty normal approximation), so it is
     only ever reached for p <= 1000, where x > 1e8 does make P = 1. special
     removed the PLIMIT branch for accuracy (its source comments on this)
     but kept the shortcut, so for p above about 1e8 every x > 1e8 now gives
     1. Removing the shortcut altogether is not enough: the continued
     fraction then overflows to NaN at x = 1e300 for every shape (16 new
     failures on `gamma_cdf`). Keeping it only for
     x > max(1e8, 2p) (where P is 1 to double precision for p >= 1e8: x is at
     least sqrt(p) = 1e4 standard deviations above the mean) removes all
     21 points; 14 points of `gammainc_large` still fail, those of group 2
     at a = 1e8 and the bulk at a = 1e9, which now reach the prefix.
   - **The log prefix** (group 2): p ln x - x - ln Gamma(p + 1) (or
     ln Gamma(p) for the fraction) is formed from three separately rounded
     terms. They reach 1.4e7 at p = 1e6 and 1.8e9 at p = 1e8 while their sum
     is only about -k^2 / 2 at x = p + k sqrt(p); the rounding of the terms
     (half an ulp of 1.4e7 is 9.3e-10) is an absolute error in the exponent,
     which `exp` turns into a relative error of the same size. libm's
     `lgamma` itself is accurate. Forming the exponent as
     p (log1p(d) - d) - ln(2 pi p) / 2 - S(p), with d = (x - p) / p (exact
     to a rounding: x - p is exact by Sterbenz's lemma) and S Stirling's
     series for ln Gamma(p + 1) - (p + 1/2) ln p + p - ln(2 pi) / 2, for
     p >= 10 and p/2 <= x <= 3p/2 (the experiment's own range, not a
     proposed fix), removes all 23 points; worst 2.5e-13 on
     `gammainc_large`. Off the grid this group starts near a = 800
     (1.4e-12), as the scan shows; it is the same cancellation as statrs's
     [finding 012](012-statrs-gamma-cdf.md).
   - **The underflow cut-off** (group 3): the result is exp(arg) only if
     arg >= ln(f64::MIN_POSITIVE) = -708.4, and 0 otherwise. exp of an
     argument between -744.4 and -708.4 is a subnormal number, so every
     subnormal result is flushed to 0. Lowering the cut-off removes both
     points; with all three changes every point of the three datasets
     passes (worst 8.2e-13).

   No call reached the checks' 5 s time limit, but neither loop has an
   iteration limit, so the checks keep it.
