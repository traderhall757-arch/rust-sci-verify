#!/usr/bin/env python3
"""Guardrail 9: a fix tested here is exactly the patch offered upstream.

A fix for a target crate lives in two places:

- `fixes/<crate>/NNNN-*.patch`: a `git format-patch` against the pinned
  release, ready to submit upstream (only with the owner's approval);
- `crates/<crate>-fix`: the patched items copied verbatim, so the harness can
  verify them against the reference data.

The fix crate declares the link in its Cargo.toml:

    [package.metadata.rsv.fix]
    upstream = "puruspe"
    patch = "fixes/puruspe/0001-....patch"
    items = { "src/beta.rs" = ["betai", ...], "src/lib.rs" = ["EPS", ...] }

This script applies the patch to a copy of the pinned upstream source (from
`cargo metadata`), and checks that every listed item (fn, const or static)
is identical, character for character, in the patched file and in the fix
crate. With --run-tests it also runs the patched crate's own test suite,
which includes the regression test that the patch adds.

Usage: scripts/check_fixes.py [--run-tests]   (from the repository root)
"""

import json
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ITEM = r"^(?:pub(?:\([\w ]+\))? )?(?:fn|const|static) {name}\b"


def extract_item(source: str, name: str) -> str | None:
    """The text of the top-level fn/const/static `name`, from its first line
    to the closing `}` or `];` at column 0 (or the same line if it ends with
    `;`). None if absent."""
    lines = source.splitlines()
    start = re.compile(ITEM.format(name=re.escape(name)))
    for i, line in enumerate(lines):
        if not start.match(line):
            continue
        if line.rstrip().endswith(";"):
            return line
        for j in range(i + 1, len(lines)):
            if lines[j] in ("}", "];"):
                return "\n".join(lines[i : j + 1])
        return None
    return None


def fix_crates(root: Path) -> list[tuple[str, Path, dict]]:
    """(package name, crate dir, fix metadata) for every workspace fix crate."""
    members = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["members"]
    out = []
    for m in members:
        package = tomllib.loads((root / m / "Cargo.toml").read_text())["package"]
        fix = package.get("metadata", {}).get("rsv", {}).get("fix")
        if fix is not None:
            out.append((package["name"], root / m, fix))
    return out


def upstream_dirs(root: Path) -> dict[str, Path]:
    """Source directory of every registry package in the locked dependency graph."""
    meta = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--locked"],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    )
    return {
        p["name"]: Path(p["manifest_path"]).parent
        for p in json.loads(meta.stdout)["packages"]
        if (p.get("source") or "").startswith("registry+")
    }


def apply_patch(upstream: Path, patch: Path, dest: Path) -> str | None:
    """Copies `upstream` to `dest` and applies `patch`; an error message on failure."""
    shutil.copytree(upstream, dest)
    out = subprocess.run(
        ["git", "apply", "--whitespace=nowarn", str(patch)],
        cwd=dest,
        capture_output=True,
        text=True,
        check=False,
    )
    return None if out.returncode == 0 else out.stderr.strip()


def compare(patched: Path, copy_source: str, items: dict[str, list[str]], label: str) -> list[str]:
    """Problems where a listed item differs between `patched` and the copy."""
    found = []
    for rel, names in items.items():
        upstream_source = (patched / rel).read_text()
        for name in names:
            want = extract_item(upstream_source, name)
            got = extract_item(copy_source, name)
            if want is None:
                found.append(f"{label}: {name} not found in patched {rel}")
            elif got is None:
                found.append(f"{label}: {name} not found in the fix crate")
            elif got != want:
                found.append(f"{label}: {name} differs from patched {rel}")
    return found


def problems(root: Path, upstreams: dict[str, Path], run_tests: bool) -> list[str]:
    found = []
    for name, crate_dir, fix in fix_crates(root):
        upstream = upstreams.get(fix["upstream"])
        if upstream is None:
            found.append(f"{name}: upstream crate {fix['upstream']!r} is not a locked dependency")
            continue
        copy_source = "\n".join(p.read_text() for p in sorted((crate_dir / "src").rglob("*.rs")))
        with tempfile.TemporaryDirectory() as tmp:
            patched = Path(tmp) / "patched"
            error = apply_patch(upstream, root / fix["patch"], patched)
            if error is not None:
                found.append(f"{name}: {fix['patch']} does not apply: {error}")
                continue
            found += compare(patched, copy_source, fix["items"], name)
            if run_tests:
                tests = subprocess.run(
                    ["cargo", "test", "-q", "--locked"],
                    cwd=patched,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if tests.returncode != 0:
                    tail = "\n".join(
                        tests.stdout.splitlines()[-20:] + tests.stderr.splitlines()[-20:]
                    )
                    found.append(f"{name}: patched upstream tests fail:\n{tail}")
    return found


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parent.parent
    found = problems(root, upstream_dirs(root), run_tests="--run-tests" in argv)
    for p in found:
        print(p)
    print(f"fix/patch check: {len(found)} problem(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
