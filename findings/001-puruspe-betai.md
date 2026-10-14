# 001: puruspe::betai wrong for large a and b, and inaccurate everywhere

- **Crate / version:** puruspe 0.4.4
- **Severity:** high (silently wrong answers: 0 instead of 0.5, 3% off), plus a
  low-severity accuracy loss at most other points
- **Status:** approved by owner (2026-10-14)
- **Affected inputs:** 336 of 639 points of `reference/data/betainc.jsonl`
  (tolerance 1e-12); the exact list is `BETAI` in
  `crates/check-puruspe/src/expected.rs`, the report is
  [reports/puruspe-0.4.4.md](../reports/puruspe-0.4.4.md).
  1. **Wrong answers** (5 points), all with a = 1e4 and b >= 1e3:
     I_x(1e4, 1e4) is 0 or -0 for x in {0.5, 0.7, 0.9, 0.99} (true: 0.5 and 1),
     and I_0.9(1e4, 1e3) is 3.1% low.
  2. **Accuracy loss** (331 points): relative errors from 1e-12 up to 5.7e-10, at
     every parameter size on the grid, including a = b = 0.1. By max(a, b):
     0.1 to 10: up to 3e-11; 100: up to 6e-11; 1000 and 1e4: up to 5.7e-10.
     SciPy's worst error on the same grid is 1.6e-13.

## Reproduction

```rust
// cargo new repro && cd repro && cargo add puruspe@=0.4.4
fn main() {
    // I_{1/2}(a, a) = 1/2 exactly, by the symmetry I_x(a, b) = 1 - I_{1-x}(b, a).
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.5)); // 0.0, should be 0.5
    println!("{:?}", puruspe::betai(1e4, 1e4, 0.7)); // -0.0, should be 1.0 (to double precision)
    println!("{:?}", puruspe::betai(1e4, 1e3, 0.9)); // 0.000596..., should be 0.000615...
    println!("{:?}", puruspe::betai(1e3, 1e3, 0.5)); // 0.5000000000864451, should be 0.5
    println!("{:?}", puruspe::betai(0.1, 0.1, 0.5)); // 0.5000000000082913, should be 0.5
}
```

## Expected

| Input (a, b, x) | Crate | mpmath (true value) | SciPy | Rel. error |
|---|---|---|---|---|
| (1e4, 1e4, 0.5) | 0.0 | 0.5 | 0.4999999999999989 | 1.0 |
| (1e4, 1e4, 0.7) | -0.0 | 1.0 | 1.0 | 1.0 |
| (1e4, 1e3, 0.9) | 0.0005960975545106217 | 0.000615295567297501 | 0.0006152955672975032 | 3.1e-2 |
| (1e3, 100, 0.9) | 0.14747917634280627 | 0.14747917625807305 | 0.1474791762580732 | 5.7e-10 |
| (1e3, 1e3, 0.5) | 0.5000000000864451 | 0.5 | 0.4999999999999992 | 1.7e-10 |
| (0.1, 0.1, 0.5) | 0.5000000000082913 | 0.5 | 0.4999999999999999 | 1.7e-11 |
| (0.1, 0.1, 0.01) | 0.32030825038093796 | 0.3203082503756263 | 0.3203082503756264 | 1.7e-11 |

## Notes

The true values rely on the identity I_{1/2}(a, a) = 1/2 where it applies, and
otherwise on the generator's positive-term series (DLMF 8.17.8), which is
tested against `mpmath.betainc` and the exact binomial sum DLMF 8.17.5
(`reference/tests/test_functions.py`). SciPy agrees with mpmath at every point.

Root cause: not yet diagnosed (milestone 3). Leads from reading the source:

- For a, b > `SWITCH` (3000) the function uses `betaiapprox` (quadrature); the
  four wrong answers at a = b = 1e4 lie in that branch.
- The prefactor is `exp(ln_gamma(a+b) - ln_gamma(a) - ln_gamma(b) + ...)`, and
  `ln_gamma` is a short Lanczos fit (g = 5, 7 terms). Its error would produce the
  roughly constant 1e-11 floor seen even at a = b = 0.1; the same floor appears
  in `gammp` ([002](002-puruspe-gammp.md)).
- The continued-fraction branch switch is `x < (a + 1)/(a + b*2)`; Numerical
  Recipes, which the code otherwise follows, uses `(a + 1)/(a + b + 2)`. The
  fifth wrong answer, (1e4, 1e3, 0.9), is affected: the written switch is 0.833,
  so it takes the `1 - ...` branch, while the Numerical Recipes switch (0.909)
  would evaluate it directly. Unverified whether this is the cause.
