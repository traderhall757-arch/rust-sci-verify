# rust-sci-verify: rules for Claude

Start every session by reading `docs/status.md`. Read `docs/roadmap.md` only when planning.

## Rules
1. Small PRs, one function family each: code, reference tests, docs together. The PR
   description explains in plain language what changed, why, how it was verified, and
   shows coverage and the mutation score.
2. Every PR has a **Checks changed** section listing any change to tolerances,
   thresholds, expected failures, conventions, tags or workflow files, with the reason.
   Write "none" if there are none. No cheap greens: never loosen a check to pass.
3. Every expected failure cites an existing `findings/NNN-crate-function.md` (CI checks).
4. Tags (pole, overflow, ...) are assigned by rules in `reference/`, never by hand.
   A finite crate value at a pole or overflow point is "crate wrong".
5. mpmath is the arbiter. A SciPy/crate disagreement is not a bug until mpmath says so.
6. Nothing speculative: code exists only if wired in, documented and reference-tested.
7. Every algorithm cites its source and license; every function documents its verified
   accuracy. Port only from clearly licensed sources (Boost, papers); not Cephes.
8. Never contact other projects (issues, PRs, comments) without the owner's explicit
   approval per item. Never publish to crates.io without approval.
9. The weekly job only reports (issue in this repo); it never changes committed data.
10. Push at the end of each meaningful step. One milestone per session; then update
    `docs/status.md`, push, and tell the owner to start a fresh session.
11. The owner merges PRs. Never commit to `main` directly.

## Commands (keep output quiet)
```sh
cargo fmt --all --check
cargo clippy -q --all-targets --locked -- -D warnings
cargo test -q --locked                         # also rewrites reports/*.md
cd reference && uv sync -q --locked && uv run -q python -m rsv_reference   # regenerate data
cd reference && uv run -q pytest -q && uv run -q ruff check . && uv run -q ruff format --check .
```
Never print whole data files or CI logs; query data with scripts, tail/grep logs.
Mutation testing runs in CI; locally only on the single function being fixed.

## Layout
- `crates/harness`: shared library (data loading, classification, reports).
- `crates/check-<crate>`: one crate per verified target crate (pins its version).
- `reference/`: pinned Python generator; data in `reference/data/<function>.jsonl`.
- `reports/`, `findings/`, `docs/` (status, roadmap, conventions, data-format).
