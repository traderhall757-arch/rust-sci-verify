# 036: statrs Chi inverse_cdf wrong for k = 1 below p = 1e-157, and at subnormal p

- **Crate / version:** statrs 0.19.1 (`statrs::distribution::Chi`, method
  `inverse_cdf` of `ContinuousCDF`; the search is
  `distribution::internal::newton_raphson_quantile`)
- **Severity:** high (silently wrong answers: `Chi::new(1).inverse_cdf(1e-300)`
  is 6.67e-162 instead of 1.25e-300, and `Chi::new(300).inverse_cdf(5e-324)`
  is 10% off)
- **Status:** approved (Claude, checklist complete, 2026-10-16)
- **Affected inputs:** 14 of 351 points of `reference/data/chi_ppf.jsonl`
  (tolerance 2e-13). The exact list is `CHI_PPF` in
  `crates/check-statrs/src/expected.rs`, the report is
  [reports/statrs-0.19.1.md](../reports/statrs-0.19.1.md). Two groups, by
  cause (section "Root cause"):
  1. **k = 1, p below about 1e-157** (2 grid points, p = 1e-300 and
     p = 5e-324; scan: every
     normal p up to 3.3e-158 fails, from 3.7e-156 up none does): for every p
     below about 6e-162 the answer is 6.668276248455232e-162 (p below about
     5e-176) or 7.202561545265052e-162, whatever p is; above that, up to
     2.5e-7 off (1.2533141155747013e-158 at p = 1e-158).
  2. **Subnormal p** (12 grid points: p = 5e-324 at every k >= 2; scan: 251
     of 298 points with p below 2.2e-308, at every k from 1 to 9988; at
     k = 1 group 1's cause gives the answer): the
     answer is far too large, by a factor of 3.0e7 at k = 2, 90,000 at
     k = 3, 24 at k = 10, and 33% (k = 100), 10% (k = 300) and 0.58%
     (k = 1e4) off.

The bulk and the upper tail are right: the scan found no failure at any
p >= 2.2e-308 except group 1. Unlike `InverseGamma::inverse_cdf`
([finding 025](025-statrs-inverse-gamma-inverse-cdf.md), the same search),
the density does not overflow and the upper-tail bracket does not miss on
this grid or in the scan. The cdf itself is
[finding 035](035-statrs-chi-cdf.md).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add statrs@=0.19.1
use statrs::distribution::{Chi, ContinuousCDF};

fn main() {
    let chi = |k| Chi::new(k).unwrap();
    println!("{:?}", chi(1).inverse_cdf(1e-300)); // should be 1.2533141373155002e-300
    println!("{:?}", chi(1).inverse_cdf(1e-200)); // should be 1.2533141373155003e-200
    println!("{:?}", chi(1).inverse_cdf(1e-158)); // should be 1.2533141373155003e-158
    println!("{:?}", chi(2).inverse_cdf(5e-324)); // should be 3.1434555694052576e-162
    println!("{:?}", chi(10).inverse_cdf(5e-324)); // should be 1.0661334955280734e-32
    println!("{:?}", chi(300).inverse_cdf(5e-324)); // should be 0.8897194915468506
    println!("{:?}", chi(10000).inverse_cdf(5e-324)); // should be 74.14756918705707
    println!("{:?}", chi(3).inverse_cdf(0.5)); // fine: should be 1.5381722544550522
}
```

This is the program that was run (release mode); it prints
`6.668276248455232e-162` (twice), `1.2533141155747013e-158`,
`9.322925914000263e-155`, `2.662405555120907e-31`, `0.9810831882742447`,
`74.57970813264653` and `1.538172254455054`.

## Expected

| k | p | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|
| 1 | 1e-300 | 6.668276248455232e-162 | 1.2533141373155002e-300 | 0.0 | 5.3e138 |
| 1 | 5e-324 | 6.668276248455232e-162 | 5e-324 | 0.0 | 3.0e146 (vs 2.2e-308) |
| 2 | 5e-324 | 9.322925914000263e-155 | 3.1434555694052576e-162 | 3.1434555694052576e-162 | 3.0e7 |
| 10 | 5e-324 | 2.662405555120907e-31 | 1.0661334955280734e-32 | 1.0661334955280863e-32 | 24 |
| 100 | 5e-324 | 0.00485351561513895 | 0.0036499223188811887 | 0.0036499223188811917 | 0.33 |
| 300 | 5e-324 | 0.9810831882742447 | 0.8897194915468506 | 0.8897194904512208 | 0.10 |
| 1e4 | 5e-324 | 74.57970813264653 | 74.14756918705707 | 74.14757034537934 | 5.8e-3 |

Second opinion: SciPy (`scipy.stats.chi.ppf`) agrees with mpmath within
the tolerance at k = 1, p = 5e-324 (0, an error of 5e-324 relative to
2.2e-308) and at every k below 300; it is 0 at
k = 1, p = 1e-300 and 1.2e-9 to 1.6e-8 off at p = 5e-324 for k >= 300.
There the closed forms are a second route: for k = 1 the quantile is
sqrt(2) erfinv(p) (1.2533141373155002e-300 = p sqrt(pi / 2) to all digits),
for k = 2 it is sqrt(-2 log(1 - p)); both are tested against the dataset's
route, and the k = 300 and 1e4 values by F(ppf(p)) = p with
`mpmath.gammainc`, in `reference/tests/test_functions.py`.

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `statrs = "=0.19.1"`
   (the program above, release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The route, sqrt(2 t) with t the gamma quantile
   at shape k/2 (`_mp_chi_ppf` on `_mp_gamma_ppf`), is tested against
   F(ppf(p)) = p (or 1 - F = 1 - p above 1/2) with `mpmath.gammainc` at 10
   points including every kind of failing one, and against the closed forms
   for k = 1 and 2, in `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 3,000 random points (Python `random.Random(20261016)`,
   shared with finding 035's scan): k = round(10^u) with u uniform in
   [0, 4]; p = 10^-v with v uniform in [308, 323.3] (10%), [0.302, 300]
   (50%), or p = 1 - 10^-v with v uniform in [0.302, 15] (40%); reference by
   the dataset's route at 30 digits. 277 fail: the 251 with subnormal p
   listed above (smallest error 2.1e-4, at k = 3094) and 26 with k = 1 and
   normal p from 4.3e-296 to 3.3e-158 (error 6.2e-10 at p = 3.3e-158, and a
   factor of 1e15 or more for p below 5e-177). None of the 2,702 other
   points fails.
5. **Root cause** (confirmed with a copy of `newton_raphson_quantile`, bit
   for bit equal to statrs at all 3,000 scan points and the 325 grid points
   with 0 < p < 1, with a switchable cdf):
   - **The search runs on `Chi::cdf`, which gives 0 well above its true
     value.** The bracket is built by halving `low` from 1 while
     `cdf(low) > p`, and the Newton step is `(cdf(x) - p) / pdf(x)`. For
     k = 1, `cdf` is 0 for x below 7.2e-162 (finding 035, causes 1 and 2;
     the halving stops at 4.4e-162, where `x * x / 2` is 1e-323, before it
     reaches a panicking x) and loses digits up to about 2e-156 (cause 3).
     So for any p below about 6e-162 the bracket is [4.4e-162, 8.9e-162]
     and the residual is -p at every x below 7.2e-162. Below about 5e-176
     the step p / 0.8 is under the relative accuracy 1e-14, and the search
     returns the midpoint 6.67e-162; above, it steps or bisects up to
     7.2e-162, the first x where `cdf` is not 0, and returns that. For k >= 2, `cdf` returns 0 wherever the true value is below
     about 1e-308 (finding 012, cause 2: `gamma_lr` returns 0 when
     `ax < -709.78`), so for subnormal p the halving and the search stop at
     the x where the true cdf is about 1e-309, not 5e-324.
   - With a cdf that does not flush (the lower series of P(k/2, x^2/2) scaled
     in log space, and x sqrt(2 / pi) for k = 1 where x^2/2 < 1e-300), the
     same search passes every point with p >= 2.2e-308, on the scan and on
     the grid (and k = 1 at p = 5e-324). What is left are 177 scan points
     and 12 grid points, all
     with subnormal p, up to 4.4% and 6.1% off: finding 025's cause 3 (a
     subnormal cdf has too few bits to locate the root; only a search on
     log F avoids it).

## Notes

`Chi::pdf` forms `x.powf(k - 1)` and `gamma(k / 2)` directly for k <= 160,
but on this grid and in the scan they stay finite wherever the search
evaluates them, so `InverseGamma`'s overflow (finding 025, cause 1) does
not occur here. SciPy computes this quantile as
`np.sqrt(2*sc.gammaincinv(.5*df, q))` (`scipy/stats/_continuous_distns.py`).
