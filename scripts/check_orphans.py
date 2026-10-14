#!/usr/bin/env python3
"""Guardrail 1: fail if any .rs file under crates/ is not compiled.

The list of compiled files comes from the compiler itself: every artifact
cargo reports for a workspace package has a dep-info file (target/.../deps/
<name>-<hash>.d) listing the sources rustc read. Only artifacts from this run
are used, so stale dep-info from earlier builds cannot hide an orphan.

Usage: scripts/check_orphans.py   (from the repository root)
"""

import json
import subprocess
import sys
from pathlib import Path

BUILDS = (
    ["cargo", "build", "-q", "--locked", "--workspace", "--all-targets"],
    ["cargo", "test", "-q", "--locked", "--workspace", "--all-targets", "--no-run"],
)


def dep_info_path(artifact: Path) -> Path:
    """target/debug/deps/libfoo-HASH.rlib -> target/debug/deps/foo-HASH.d"""
    stem = artifact.name.split(".")[0]
    stem = stem.removeprefix("lib") if artifact.suffix in {".rlib", ".rmeta", ".so"} else stem
    return artifact.parent / f"{stem}.d"


def parse_dep_info(text: str, root: Path) -> set[Path]:
    """Source paths named on the left of 'path:' lines in a make-style .d file."""
    sources = set()
    for line in text.splitlines():
        if line.endswith(":") and line.strip().endswith(".rs:"):
            p = Path(line[:-1].strip())
            sources.add((p if p.is_absolute() else root / p).resolve())
    return sources


def compiled_sources(root: Path) -> set[Path]:
    sources: set[Path] = set()
    for cmd in BUILDS:
        out = subprocess.run(
            [*cmd, "--message-format=json"], cwd=root, check=True, capture_output=True, text=True
        ).stdout
        for line in out.splitlines():
            msg = json.loads(line)
            if msg.get("reason") != "compiler-artifact" or "path+file://" not in msg["package_id"]:
                continue
            for name in msg["filenames"]:
                d = dep_info_path(Path(name))
                if d.is_file():
                    sources |= parse_dep_info(d.read_text(), root)
    return sources


def orphans(root: Path, compiled: set[Path]) -> list[Path]:
    present = {p.resolve() for p in (root / "crates").rglob("*.rs") if "target" not in p.parts}
    return sorted(present - compiled)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    found = orphans(root, compiled_sources(root))
    for p in found:
        print(f"orphan (not compiled): {p.relative_to(root)}")
    print(f"orphan check: {len(found)} orphan file(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
