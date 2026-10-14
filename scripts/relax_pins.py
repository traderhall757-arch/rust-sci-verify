#!/usr/bin/env python3
"""Weekly job only: let check crates resolve their target crate's latest release.

Rewrites `name = "=x.y.z"` to `name = "*"` in every check crate's manifest.
Runs in the CI working copy, which the weekly job never commits.

Usage: scripts/relax_pins.py   (from the repository root)
"""

import re
import sys
import tomllib
from pathlib import Path

PIN = re.compile(r'^(\s*[\w-]+\s*=\s*)"=\d+\.\d+\.\d+"', re.MULTILINE)


def relax(text: str) -> str:
    return PIN.sub(r'\1"*"', text)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    for m in tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["members"]:
        manifest = root / m / "Cargo.toml"
        text = manifest.read_text()
        if tomllib.loads(text)["package"].get("metadata", {}).get("rsv", {}).get("role") == "check":
            manifest.write_text(relax(text))
            print(f"relaxed pins in {m}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
