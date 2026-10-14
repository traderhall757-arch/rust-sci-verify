#!/usr/bin/env python3
"""Guardrail 6: every public function has a reference test.

Every workspace crate declares its role in Cargo.toml:

    [package.metadata.rsv]
    role = "harness" | "check" | "implementation"

- implementation: scientific code we write (Phases B and C). Every `pub fn`
  must be named by a reference check, i.e. appear as
  `function: "<crate>::<...>::<name>"` in some check crate's sources.
- check: verifies a target crate; must define no `pub fn` at all.
- harness: test tooling with no scientific functions; covered instead by the
  coverage and mutation-testing guardrails.

Usage: scripts/check_api_tests.py   (from the repository root)
"""

import re
import sys
import tomllib
from pathlib import Path

ROLES = {"harness", "check", "implementation"}
# `pub fn` only: pub(crate) items are internal, and the unreachable_pub lint
# guarantees that every plain `pub` item is part of the public API.
PUB_FN = re.compile(r"^\s*pub\s+(?:const\s+|async\s+|unsafe\s+)*fn\s+(\w+)", re.MULTILINE)
CHECKED = re.compile(r'function:\s*"([\w:]+)"')


def crates(root: Path) -> list[tuple[str, str, Path]]:
    """(package name, role, crate dir) for every workspace member."""
    members = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["members"]
    out = []
    for m in members:
        manifest = tomllib.loads((root / m / "Cargo.toml").read_text())
        role = manifest.get("package", {}).get("metadata", {}).get("rsv", {}).get("role")
        out.append((manifest["package"]["name"], role, root / m))
    return out


def sources(crate_dir: Path) -> str:
    return "\n".join(p.read_text() for p in sorted((crate_dir / "src").rglob("*.rs")))


def problems(root: Path) -> list[str]:
    found = []
    all_crates = crates(root)
    checked = set()
    for _, role, d in all_crates:
        if role == "check":
            checked |= set(CHECKED.findall(sources(d)))
    for name, role, d in all_crates:
        if role not in ROLES:
            found.append(f"{name}: [package.metadata.rsv] role must be one of {sorted(ROLES)}")
            continue
        fns = PUB_FN.findall(sources(d))
        if role == "check" and fns:
            found.append(f"{name}: check crates must not define pub fn ({', '.join(fns)})")
        if role == "implementation":
            lib = name.replace("-", "_")
            for f in fns:
                if not any(c.startswith(f"{lib}::") and c.endswith(f"::{f}") for c in checked):
                    found.append(f"{name}: pub fn {f} has no reference check")
    return found


def main() -> int:
    found = problems(Path(__file__).resolve().parent.parent)
    for p in found:
        print(p)
    print(f"api/reference-test check: {len(found)} problem(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
