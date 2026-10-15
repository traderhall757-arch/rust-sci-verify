# Reference data format (schema 1)

One file per mathematical function: `reference/data/<function>.jsonl`, written by
`python -m rsv_reference` and read by `crates/harness`. Files are committed so
`cargo test` needs no Python. CI regenerates them and fails if they differ
(on pull requests only the datasets they can have changed; see `docs/guardrails.md`).

## Layout: JSON Lines

The **first line** is the header:

```json
{"schema":1,"function":"erfcx","description":"...","inputs":["x"],"outputs":["y"],
 "points":361,"tolerance":{"rel":1e-13,"reason":"..."},"grid":"...",
 "mpmath":{"version":"1.4.1","dps":[50,100]},
 "scipy":{"version":"1.18.1","function":"scipy.special.erfcx"},
 "numpy":"2.5.4","python":"3.13"}
```

Every **following line** is one point:

```json
{"in":{"x":[-30.0]},"mpmath":{"y":["inf"]},"scipy":{"y":["inf"]},"tags":["overflow"]}
```

- `in`, `mpmath`, `scipy` map each input or output name to a **list** of numbers.
  Scalars are lists of length 1. Phase B (spectral analysis) uses the same layout
  with longer lists, e.g. `"in":{"x":[...4096 samples...],"fs":[1000.0]}` and
  outputs `"f"` and `"Pxx"`. The harness compares lists element by element and
  reports the worst element.
- Numbers are written as the **shortest decimal that reads back to the same
  double** (Python `repr`). The Rust reader parses with `serde_json`'s
  `float_roundtrip` feature, so every value is bit-exact.
- Non-finite values are strings: `"nan"`, `"inf"`, `"-inf"`.
- `tags` is omitted when empty. See [conventions.md](conventions.md).
- mpmath values are the 50-digit result **correctly rounded** to a double
  (round to nearest, ties to even, subnormals included).

## Determinism

Grids are built from exact decimal strings or exact binary fractions, never from
platform `pow()`, and there is no randomness. The SciPy column can differ in
the last bit between platforms, so CI (`ubuntu-24.04`, pinned Python
environment from `reference/uv.lock`) is the source of truth.

## Adding a function

Add a `Spec` to `reference/rsv_reference/functions.py` (grid, mpmath and SciPy
callables, tolerance with its reason, tag rules), test its mpmath route against
an independent one in `reference/tests/test_functions.py`, regenerate, then add
a `Check` to each relevant `crates/check-<crate>/src/verify.rs`. Multi-input
functions list their inputs in the header (e.g. `"inputs":["a","b","x"]`); the
harness flattens them in that order to match expected failures.
