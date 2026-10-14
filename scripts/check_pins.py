#!/usr/bin/env python3
"""Guardrail 8: every version that affects results is pinned exactly.

- rust-toolchain.toml names an exact release (x.y.z);
- Cargo.lock and reference/uv.lock are committed;
- reference/pyproject.toml pins every dependency with ==;
- check crates pin their target crate with =x.y.z (new releases are tried by
  the weekly job, then adopted in a pull request);
- CI tools installed with taiki-e/install-action carry an @version.

Usage: scripts/check_pins.py   (from the repository root)
"""

import re
import subprocess
import sys
import tomllib
from pathlib import Path

EXACT = re.compile(r"^\d+\.\d+\.\d+$")


def tracked(root: Path, path: str) -> bool:
    out = subprocess.run(
        ["git", "ls-files", path], cwd=root, capture_output=True, text=True, check=False
    )
    return out.stdout.strip() == path


def problems(root: Path) -> list[str]:
    found = []
    channel = tomllib.loads((root / "rust-toolchain.toml").read_text())["toolchain"]["channel"]
    if not EXACT.match(channel):
        found.append(f"rust-toolchain.toml: channel {channel!r} is not an exact x.y.z release")
    for lock in ("Cargo.lock", "reference/uv.lock"):
        if not tracked(root, lock):
            found.append(f"{lock}: not committed")
    project = tomllib.loads((root / "reference/pyproject.toml").read_text())
    groups = project.get("dependency-groups", {}).values()
    for dep in [*project["project"]["dependencies"], *[d for g in groups for d in g]]:
        if not re.match(r"^[\w.-]+==[\w.]+$", dep):
            found.append(f"reference/pyproject.toml: {dep!r} is not pinned with ==")
    members = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["members"]
    for m in members:
        manifest = tomllib.loads((root / m / "Cargo.toml").read_text())
        if manifest["package"].get("metadata", {}).get("rsv", {}).get("role") != "check":
            continue
        deps = {**manifest.get("dependencies", {}), **manifest.get("dev-dependencies", {})}
        for name, spec in deps.items():
            if isinstance(spec, dict) and spec.get("workspace"):
                continue
            version = spec if isinstance(spec, str) else spec.get("version", "")
            if not (version.startswith("=") and EXACT.match(version[1:])):
                found.append(f"{m}: {name} = {version!r} is not pinned with =x.y.z")
    for wf in sorted((root / ".github/workflows").glob("*.yml")):
        for tool in re.findall(r"^\s*tool:\s*(\S+)", wf.read_text(), re.MULTILINE):
            for t in tool.split(","):
                if "@" not in t:
                    found.append(f"{wf.name}: install-action tool {t!r} has no @version")
    return found


def main() -> int:
    found = problems(Path(__file__).resolve().parent.parent)
    for p in found:
        print(p)
    print(f"pin check: {len(found)} problem(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
