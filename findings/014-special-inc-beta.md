# 014: special Beta::inc_beta loses up to 8 digits for large parameters

- **Crate / version:** special 0.14.2 (`special::Beta::inc_beta`, the
  regularized incomplete beta function I_x(a, b); default features)
- **Severity:** low (tiny inaccuracy: up to 4.0e-9 relative on the grid,
  6.0e-10 off it; no wrong answers, NaN or hangs)
- **Status:** approved (Claude, checklist complete, 2026-10-14)
- **Affected inputs:** 182 of 969 points (tolerance 1e-12): 29 of 639 in
  `reference/data/betainc.jsonl`, 96 of 144 in `betainc_large.jsonl`, 57 of
  186 in `betainc_switch.jsonl`. The exact lists are `INC_BETA`,
  `INC_BETA_LARGE` and `INC_BETA_SWITCH` in
  `crates/check-special/src/expected.rs`; the report is
  [reports/special-0.14.2.md](../reports/special-0.14.2.md). Every failing
  point has max(a, b) >= 100 on the grid (>= 283.8 off it). Errors grow with
  the parameters: up to 6.6e-11 at a or b = 1e4 (`betainc`), 1.8e-10 at
  a = b = 1e5 (`betainc_large`), and 4.0e-9 for a = 1e5 with small b and x
  close to 1 (`betainc_switch`).

`inc_beta` is called as its documentation asks, with `ln_beta` precomputed
by `a.ln_beta(b)`. Its inverse has its own, worse problems
([finding 015](015-special-inv-inc-beta.md)).

## Reproduction

```rust
// cargo new repro && cd repro && cargo add special@=0.14.2
use special::Beta;

fn inc_beta(a: f64, b: f64, x: f64) -> f64 {
    x.inc_beta(a, b, a.ln_beta(b)) // ln_beta precomputed, as the documentation asks
}

fn main() {
    println!("{:?}", inc_beta(1e5, 0.01, 0.99999999)); // should be 0.06142939379322011
    println!("{:?}", inc_beta(1e5, 1e5, 0.5)); // should be 0.5 (symmetry)
    println!("{:?}", inc_beta(1e4, 1e4, 0.5)); // should be 0.5 (symmetry)
    println!("{:?}", inc_beta(1e4, 0.1, 0.9999900000999989)); // should be 0.17245201791471385
    println!("{:?}", inc_beta(1e4, 5.0, 0.99)); // should be 9.76189855487456e-38
    let t = std::time::Instant::now();
    let v = inc_beta(1e5, 0.05, 0.999999); // should be 0.08874255079357986
    println!("{v:?} in {:?}", t.elapsed());
}
```

This is the program that was run (release mode). It prints the "Crate"
column below; the last call took 127 ms.

## Expected

| a | b | x | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|---|---|
| 1e5 | 0.01 | 0.99999999 | 0.06142939404049452 | 0.06142939379322011 | 0.06142939379322009 | 4.0e-9 |
| 1e5 | 0.05 | 0.999999 | 0.08874255086682581 | 0.08874255079357986 | 0.08874255079357986 | 8.3e-10 |
| 1e5 | 1e5 | 0.5 | 0.5000000000909487 | 0.5 | 0.49999999999999895 | 1.8e-10 |
| 1e4 | 0.1 | 0.9999900000999989 | 0.17245201792601522 | 0.17245201791471385 | 0.17245201791471385 | 6.6e-11 |
| 1e4 | 5 | 0.99 | 9.761898554583739e-38 | 9.76189855487456e-38 | 9.76189855487456e-38 | 3.0e-11 |
| 1e4 | 1e4 | 0.5 | 0.4999999999859591 | 0.5 | 0.4999999999999989 | 2.8e-11 |

Second opinion: SciPy (`scipy.special.betainc`) agrees with mpmath within
2.2e-15 at every row; at a = b, x = 1/2 the true value is exactly 1/2 by the
symmetry I_x(a, b) = 1 - I_(1-x)(b, a).

## Confirmation checklist

1. **Reproduced** in a fresh `cargo new` project with `special = "=0.14.2"`
   (default features, as in the check crate: `no_std`, which takes `ln`,
   `exp` and `lgamma` from the `libm` crate; the program above, built in
   release mode).
2. **Arbiter checked**: the generator requires mpmath at 50 and 100 digits to
   round to the same double. The three datasets use the positive-term series
   DLMF 8.17.8, tested against `mpmath.betainc`, the exact binomial sum and
   the continued fraction DLMF 8.17.22 in
   `reference/tests/test_functions.py`.
3. **Second opinion**: see "Expected".
4. **Off-grid**: 10,000 random points (Python `random.Random(20261014)`):
   a, b = 10^u with u uniform in [-1, 5]; x uniform in (0, 1) for half the
   points and x = SciPy's `betaincinv(a, b, p)` with p uniform in (0, 1) for
   the other half (so that I is not extreme); each rounded to 6 significant
   digits. The crate was compared with SciPy, and mpmath (the continued
   fraction at 50 and 100 digits; they agreed at every point) decided
   wherever the two differed by more than 1e-13. 2,362 points fail
   (> 1e-12), all with max(a, b) >= 283.8: 2,305 of the 5,536 points with
   max(a, b) >= 1000 and 1,486 of the 3,021 with max(a, b) >= 1e4. Worst
   6.0e-10 (a = 68082.5, b = 64755.7, x = 0.509026); below max(a, b) = 1000
   at most 3.5e-12. No call failed to return within 20 s.
5. **Root cause** (confirmed with a copy of `inc_beta`, scratch, not
   committed, that returns the crate's result bit for bit at all 969 grid
   points and can swap parts for exact ones). `inc_beta` is AS 63 (Soper's
   reduction series) with the AS R19 / AS 109 prefactor:
   - **Prefactor** (139 of the 182 points): the series sum is multiplied by
     `exp(p ln x + (q - 1) ln(1 - x) - ln_beta) / p`. The three terms of the
     exponent are each rounded to a double, and `ln_beta` =
     ln Gamma(a) + ln Gamma(b) - ln Gamma(a + b) is itself a difference of
     such rounded terms (libm's `lgamma`), up to 2.2e6 in magnitude
     (ln Gamma(2e5)), while the exponent's value is small. So the exponent
     has an absolute error of up to about 2.2e6 * 2^-53 = 2.5e-10, which
     `exp` turns into a relative error of the same size. Where the series is evaluated for the complement
     (AS 63 swaps the parameters when a < (a + b) x) and the result is
     1 - I', the error is magnified by I' / (1 - I'), up to about 15 at
     a = 1e5, b = 0.01. With `ln_beta` exact (mpmath, rounded) the failures
     drop from 182 to 64; with the whole prefactor exact, to 43 (worst
     9.0e-10).
   - **Series cut-off** (12 of the remaining 43): the series stops when a
     term is below 1e-15 times the sum, but for small b with x near 1 the
     terms shrink so slowly that the omitted tail is much larger than the
     last term: at a = 1e5, b = 0.05, x = 0.999999 the series runs for
     14.8 million terms (hence the 127 ms above) and still stops early.
     With the prefactor exact and the cut-off at 1e-30, 31 points remain,
     worst 5.1e-11.
   - **Lead, not confirmed** (the last 31, up to 1.1e-11 at a = b = 1e5 near
     x = 1/2 and 5.1e-11 for small b): rounding errors accumulated over
     the long series (277 to 47 million terms at these points).
