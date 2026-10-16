# 012: statrs Gamma and ChiSquared cdf and sf lose digits for large shape

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Gamma` and
  `statrs::distribution::ChiSquared`, methods `cdf` and `sf` of
  `ContinuousCDF`; both call `statrs::function::gamma::gamma_lr` /
  `gamma_ur`; `statrs::distribution::Erlang` calls `Gamma`'s, and
  `statrs::distribution::InverseGamma` calls them at 1/x, and
  `statrs::distribution::Chi` at x^2/2)
- **Severity:** low (an accuracy loss: up to 1.2e-9 on the grid and
  3.4e-9 off it, about 9 correct digits), plus subnormal results
  flushed to 0
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 22 of 538 points of `reference/data/gamma_cdf.jsonl`
  and 21 of 538 points of `reference/data/gamma_sf.jsonl` (tolerance 1e-12,
  the `gammainc` tolerance); the exact lists are `GAMMA_CDF` and `GAMMA_SF`
  in `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Only rate 1 is
  checked. `ChiSquared::new(k)` is `Gamma::new(k / 2, 0.5)` and multiplies x
  by 0.5 exactly, so `ChiSquared` with df = 2a at 2x fails at exactly the
  same points with the same values (checked on both datasets). `Erlang`
  (whole-number shape) calls `Gamma`'s `cdf` and `sf` directly; every point
  of both lists has a whole-number shape, so `Erlang` fails at all of them
  (checked on the 409 whole-number-shape points of each dataset, added
  2026-10-15). `InverseGamma::cdf` is `gamma_ur(a, 1/x)` and `sf` is
  `gamma_lr(a, 1/x)`, so it fails at 8 of 531 points of
  `reference/data/invgamma_cdf.jsonl` and 10 of 531 of
  `reference/data/invgamma_sf.jsonl` for the same two causes (lists
  `INVGAMMA_CDF` and `INVGAMMA_SF`; added 2026-10-15): 16 points with
  a >= 1000, worst 2.5e-11, and the subnormal results at a = 5, x = 1/750
  (`cdf`, true 2.5e-316) and a = 50, x = 1/1e-5 (`sf`, true 3.3e-315). Rounding
  1/x costs at most 3.4e-13 on that grid (its tolerance reason). Its
  panics for tiny x are finding 024. `Chi::cdf` is `gamma_lr(k/2, x^2/2)` and
  `sf` is `gamma_ur(k/2, x^2/2)`, so it fails at 6 of 514 points of
  `reference/data/chi_cdf.jsonl` and 7 of 514 of
  `reference/data/chi_sf.jsonl` for the same two causes (lists `CHI_CDF` and
  `CHI_SF`; added 2026-10-16): 9 points at k = 1e4 (shape 5e3), worst
  8.6e-12, and subnormal results returned as 0 (`cdf` at k = 100, true
  3.3e-315; `sf` at k = 10, 160 and 161, true 2.5e-316 to 2.2e-314).
  Rounding x * x costs at most 2.4e-13 on that grid (its tolerance reason);
  in finding 035's scan the errors start at k = 1454. Its panics where
  x * x under- or overflows are finding 035. Called directly (added 2026-10-15),
  `statrs::function::gamma::gamma_lr` fails at 1 of 192 points of
  `reference/data/gammainc.jsonl` (a = 200, x = 2: true 2.8e-316, returned
  as 0; list `GAMMA_LR`; its panics at x = 0 are finding 032) and at 14 of
  51 points of `reference/data/gammainc_large.jsonl` (7 each at a = 1e8
  and 1e9, across the bulk and both tails down to 2e-198; list
  `GAMMA_LR_LARGE`). There cause 1 reaches 2.6e-7 (a = 1e8) and 4.1e-6
  (a = 1e9): the spacing of doubles near a ln a (1.8e9 and 2.1e10) is 2.4e-7
  and 3.8e-6, the absolute error that rounding the terms of `ax` leaves.
  Two groups,
  by cause (section "Root cause"):
  1. **Large shape** (41 points, all with a >= 1000): the error grows like
     the shape, worst per shape 1.4e-12 (a = 1000), 2.4e-11 (1e4), 4.3e-10
     (1e5) and 1.2e-9 (1e6), across the bulk and both tails.
  2. **Subnormal results flushed to 0** (2 points): `cdf` at a = 50,
     x = 1e-5 (true 3.3e-315) and `sf` at a = 5, x = 750 (true 2.5e-316).
     Every result below about 1e-308 is returned as 0 (or its complement 1).

The quantile `inverse_cdf` has separate, more severe problems:
[finding 013](013-statrs-gamma-inverse-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{ChiSquared, ContinuousCDF, Gamma};

fn main() {
    let g = |a| Gamma::new(a, 1.0).unwrap();
    println!("{:?}", g(1e6).cdf(997000.0)); // should be 0.0013381041673135997
    println!("{:?}", g(1e6).sf(1e6)); // should be 0.4998670192391274
    println!("{:?}", g(1e5).cdf(1e5)); // should be 0.5004205221103651
    println!("{:?}", g(1e4).cdf(9000.0)); // should be 2.073299202433928e-25
    println!("{:?}", g(1e4).sf(13000.0)); // should be 4.712347177485643e-166
    println!("{:?}", g(1e3).cdf(905.1316701949486)); // should be 0.0010012389950642836
    println!("{:?}", g(50.0).cdf(1e-5)); // 0, should be 3.28791718e-315
    println!("{:?}", g(5.0).sf(750.0)); // 0, should be 2.52052905e-316
    // ChiSquared(k) is Gamma(k/2, rate 1/2): the same numbers.
    println!("{:?}", ChiSquared::new(2e6).unwrap().cdf(1994000.0)); // should be 0.0013381041673135997
}
```

This is the program that was run (release mode); it gives the "Crate"
column below.

## Expected

| a | x | Method | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 1e6 | 997000 | `cdf` | 0.001338104168893909 | 0.0013381041673135997 | 0.0013381041673135986 | 1.2e-9 |
| 1e6 | 1e6 | `sf` | 0.4998670195812436 | 0.4998670192391274 | 0.4998670192391274 | 6.8e-10 |
| 1e5 | 1e5 | `cdf` | 0.5004205220662236 | 0.5004205221103651 | 0.5004205221103651 | 8.8e-11 |
| 1e4 | 9000 | `cdf` | 2.0732992024016744e-25 | 2.073299202433928e-25 | 2.0732992024338907e-25 | 1.6e-11 |
| 1e4 | 13000 | `sf` | 4.7123471773705207e-166 | 4.712347177485643e-166 | 4.7123471774867427e-166 | 2.4e-11 |
| 1000 | 905.1316701949486 | `cdf` | 0.0010012389950632714 | 0.0010012389950642836 | 0.0010012389950642825 | 1.0e-12 |
| 50 | 1e-5 | `cdf` | 0.0 | 3.28791718e-315 | 0.0 | 1.5e-7 |
| 5 | 750 | `sf` | 0.0 | 2.52052905e-316 | 0.0 | 1.1e-8 |

(x = 905.1316701949486 is a - 3 sqrt(a) for a = 1000. Errors of subnormal
results are relative to 2.2e-308, per `docs/conventions.md`.)

Second opinion: SciPy (`scipy.stats.gamma`) agrees with mpmath within
2.3e-13 at every row except the two subnormal ones, where it returns 0 as
well. There, and at every integer shape, the finite sum
Q(a, x) = e^-x sum_(k<a) x^k / k! (DLMF 8.4.10) is a second route; it is
tested against the dataset's route at a = 5, x = 750 and up to a = 1e6
(`reference/tests/test_functions.py`). For a = 50, x = 1e-5 the first two terms
of the series, x^50 e^-x / 50! (1 + x / 51) = 3.28791718e-315, already give
the double.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, built in release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. `mpmath.gammainc` does not converge for
   a = 1e6 near the mean, so the dataset uses the positive-term series
   DLMF 8.7.1 for x < a + 1 and Legendre's continued fraction DLMF 8.9.2
   above, with 60 extra bits. Both are tested against `mpmath.gammainc` (up
   to a = 1000), the finite sum DLMF 8.4.10 (integer a up to 1e6), the
   closed forms for a = 1/2 (erf, erfc) and a = 1 (exp), and P + Q = 1, in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 20,000 random points (Python `random.Random(20261014)`):
   a = 10^u with u uniform in [-2, 6]; x = a + z sqrt(a) with z uniform in
   [-40, 40] (60%; where that is not positive, x = a 10^v with v uniform in
   [-30, 0]) or x = a 10^v with v uniform in [-30, 4] (40%); reference by
   the dataset's route at 30 digits. `cdf` fails 1,780 (8.9%) and `sf`
   2,009 (10.0%). Of these, 1,725 and 1,995 have a normal true value, all
   with a >= 754; worst per decade of a: 1.6e-12 (a from 100 to 1000),
   2.1e-11, 2.4e-10 and 3.4e-9 (1e5 to 1e6). The other 55 and 14 are
   subnormal results returned as 0. All 3,789 pass with the two fixes
   below (worst error left 1.4e-13), so every failure has one of the causes
   below.
5. **Root cause** (confirmed by experiment): `gamma_lr` and `gamma_ur`
   (`src/function/gamma.rs`) are the Cephes-style series and continued
   fraction, scaled by `exp(ax)` with
   `ax = a * x.ln() - x - ln_gamma(a)`. A copy of both functions (scratch,
   not committed; it reproduces statrs bit for bit on all 1,044 grid
   evaluations with x > 0) with switchable fixes shows:
   - **Cause 1:** for large a the three terms of `ax` are each about
     a ln a (1.3e7 at a = 1e6) while their sum is small (about
     -(x - a)^2 / (2a) - ln(2 pi a) / 2). Each term is rounded to a double
     whose spacing is 1.9e-9 at 1.3e7, so `ax` carries an absolute error of
     that size, and `exp(ax)` turns it into the same relative error in P.
     statrs's `ln_gamma` is not the problem: it is within about one ulp of
     ln Gamma(a) for a = 1e3 to 1e6 (measured). Supplying the correctly
     rounded `ax` (from mpmath) to the copy brings every one of the 41
     points within 2.8e-13. A method that forms x - a exactly and expands
     in (x - a) / a (Temme's uniform expansion, DLMF 8.12, used by Boost
     and SciPy for large a) avoids the cancellation.
   - **Cause 2:** both functions return exactly 0 or 1 when
     `ax < -709.78` (where `exp(ax)` would be below the smallest normal
     double), although the result `exp(ax) * sum / a` is still a subnormal
     number. Scaling in log space instead, `exp(ax + ln(sum / a))`, gives
     both subnormal points exactly.

## Notes

`gamma_lr` and `gamma_ur` are also public (`statrs::function::gamma`), so
direct callers get the same errors; this project checks them only through
the distributions. Severity is low: the largest error is about 1e-9, and the
flushed values are below 1e-308. The gamma quantile (finding 013) inherits
both causes.
