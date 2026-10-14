# Contributing

This project verifies Rust scientific crates against mpmath (the arbiter) and
SciPy (the practical reference). The rules below apply to everyone, human or AI.

## Pull requests

- **Small**: one function family per PR. Code, reference tests, docs and
  examples go together.
- **Plain-language description**: what changed, why, and how it was verified,
  with coverage and the mutation score.
- **"Checks changed" section, always**: list every change to tolerances,
  thresholds, expected failures, conventions, tags or workflow files, with the
  reason. Write "none" if there are none.
- **No cheap greens**: never loosen a tolerance, add an expected failure or
  change a tag rule to make CI pass. Every expected failure cites an existing
  `findings/` file (CI enforces this).
- **Nothing speculative**: code exists only if it is used, documented and
  reference-tested.
- **Correct beats concise**: every algorithm cites its source and its license;
  every function documents its verified accuracy. Port only from clearly
  licensed sources (Boost, published papers); not from Cephes.

## Before pushing

```sh
cargo fmt --all --check
cargo clippy -q --all-targets --locked -- -D warnings
cargo test -q --locked
cd reference && uv sync --locked && uv run python -m rsv_reference && git diff --stat
python3 scripts/check_orphans.py && python3 scripts/check_api_tests.py && python3 scripts/check_pins.py
```

CI also runs coverage (at least 90% of lines), `cargo machete`, and mutation
testing on the changed code. See [docs/guardrails.md](docs/guardrails.md).

Cloud sessions of Claude Code run `.claude/hooks/session-start.sh`, which
installs the pinned Python environment and fetches Rust dependencies.

`cargo test` rewrites `reports/*.md`; commit them with your change.

## Upstream projects

We never open issues, pull requests or comments on other projects without the
owner's explicit approval for each item. Approved reports follow that project's
contributing rules, disclose AI assistance, and stay minimal and reproducible.
Nothing is published to crates.io without approval.

## License

Contributions are dual licensed under MIT OR Apache-2.0, as described in the README.
