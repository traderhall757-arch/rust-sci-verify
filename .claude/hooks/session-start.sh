#!/bin/bash
# SessionStart hook for Claude Code cloud sessions: installs the pinned Python
# reference environment and prefetches Rust dependencies so `cargo test`,
# the generator and the guardrail scripts work from the first turn.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

# Rust: the toolchain is pinned by rust-toolchain.toml; rustup installs it on
# first use. Fetch crates so later builds need no network round trips.
cargo fetch -q --locked

# Python: exact versions from reference/uv.lock (uv sync is idempotent).
if ! command -v uv >/dev/null 2>&1; then
  pip install -q uv
fi
(cd reference && uv sync -q --locked)
