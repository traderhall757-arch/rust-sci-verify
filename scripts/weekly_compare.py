#!/usr/bin/env python3
"""Weekly job: compare a run on the latest crate and SciPy releases with the
committed results, and write a markdown report of the differences.

Compared:
- reference data: per dataset, how many points' mpmath or SciPy values changed
  (from `git show HEAD:<file>` against the regenerated file);
- reports: the per-function summary rows of each crate, matched by crate name
  (versions may differ);
- the `cargo test` outcome: unexpected failures or expected failures that no
  longer fail.

Exit status 0 = no differences, 10 = differences (report written), other = error.

Usage: scripts/weekly_compare.py COMMITTED_REPORTS_DIR TEST_STATUS TEST_LOG OUT_MD
"""

import json
import re
import subprocess
import sys
from pathlib import Path

CHANGED = 10
ROW = re.compile(r"^\| `([^`]+)` \|(.*)\|$", re.MULTILINE)
REPORT = re.compile(r"^(?P<crate>.+)-(?P<version>\d+\.\d+\.\d+.*)\.md$")


def points(text: str) -> list[dict]:
    return [json.loads(line) for line in text.splitlines()[1:]]


def data_changes(root: Path) -> list[str]:
    out = []
    for path in sorted((root / "reference/data").glob("*.jsonl")):
        rel = path.relative_to(root).as_posix()
        old = subprocess.run(
            ["git", "show", f"HEAD:{rel}"], cwd=root, capture_output=True, text=True, check=False
        )
        if old.returncode != 0:
            out.append(f"- `{rel}`: new dataset")
            continue
        a, b = points(old.stdout), points(path.read_text())
        if len(a) != len(b):
            out.append(f"- `{rel}`: point count {len(a)} -> {len(b)}")
            continue
        mp = sum(x["mpmath"] != y["mpmath"] for x, y in zip(a, b, strict=True))
        sp = sum(x["scipy"] != y["scipy"] for x, y in zip(a, b, strict=True))
        if mp or sp:
            out.append(f"- `{rel}`: {mp} mpmath value(s) and {sp} SciPy value(s) changed")
    return out


def tool_versions(root: Path) -> list[str]:
    """Reference tool versions: committed header vs regenerated header."""
    files = sorted((root / "reference/data").glob("*.jsonl"))
    if not files:
        return []
    rel = files[0].relative_to(root).as_posix()
    old = subprocess.run(
        ["git", "show", f"HEAD:{rel}"], cwd=root, capture_output=True, text=True, check=False
    ).stdout
    new_h = json.loads(files[0].read_text().splitlines()[0])
    old_h = json.loads(old.splitlines()[0]) if old else {}

    def get(h: dict, tool: str) -> str:
        v = h.get(tool, "?")
        return v["version"] if isinstance(v, dict) else v

    return [
        f"- {tool}: committed {get(old_h, tool)}, latest {get(new_h, tool)}"
        for tool in ("scipy", "numpy", "mpmath")
    ]


def summary_rows(reports: Path) -> dict[str, tuple[str, dict[str, str]]]:
    """crate -> (version, {function: summary row})"""
    out = {}
    for f in sorted(reports.glob("*.md")):
        m = REPORT.match(f.name)
        if m:
            rows = {fn: rest.strip() for fn, rest in ROW.findall(f.read_text())}
            out[m["crate"]] = (m["version"], rows)
    return out


def report_changes(committed: Path, latest: Path) -> list[str]:
    old, new = summary_rows(committed), summary_rows(latest)
    out = []
    for crate, (version, rows) in sorted(new.items()):
        old_version, old_rows = old.get(crate, ("none", {}))
        changed = sorted(
            fn for fn in rows.keys() | old_rows.keys() if rows.get(fn) != old_rows.get(fn)
        )
        if not changed:
            continue
        out.append(f"### {crate}: committed {old_version}, latest {version}")
        out.append("")
        out.append(
            "| Function | Committed (points, pass, wrong, both, SciPy wrong, conv., skipped, worst, tol.) | Latest |"
        )
        out.append("|---|---|---|")
        out += [f"| `{fn}` | {old_rows.get(fn, '-')} | {rows.get(fn, '-')} |" for fn in changed]
        out.append("")
    return out


def test_problems(status: int, log: str) -> list[str]:
    if status == 0:
        return []
    lines = [ln.strip() for ln in log.splitlines() if "unexpected" in ln or "no longer fails" in ln]
    if not lines:  # e.g. a new release no longer compiles against the checks
        tail = log.splitlines()[-15:]
        return [f"`cargo test` failed (exit {status}). Last lines:", "", "```", *tail, "```"]
    shown = [f"- `{ln}`" for ln in lines[:30]]
    more = [f"- ... and {len(lines) - 30} more"] if len(lines) > 30 else []
    return [f"`cargo test` failed (exit {status}).", "", *shown, *more]


def main(argv: list[str]) -> int:
    committed, status, log, out_md = Path(argv[0]), int(argv[1]), Path(argv[2]), Path(argv[3])
    root = Path(__file__).resolve().parent.parent
    data = data_changes(root)
    reports = report_changes(committed, root / "reports")
    tests = test_problems(status, log.read_text())
    if not (data or reports or tests):
        print("weekly compare: no differences from committed results")
        return 0
    body = ["# Weekly verification: results differ from committed data", ""]
    body += [
        "Latest releases of the target crates and of SciPy/NumPy/mpmath were run in CI.",
        "Nothing was committed; adopt changes through a normal pull request.",
        "",
    ]
    body += ["## Reference tool versions", "", *tool_versions(root), ""]
    if tests:
        body += ["## Test outcome", "", *tests, ""]
    if reports:
        body += ["## Report summaries", "", *reports]
    if data:
        body += ["## Reference data", "", *data, ""]
    out_md.write_text("\n".join(body))
    print(
        f"weekly compare: differences found ({len(data)} dataset(s), {len(reports)} report line(s))"
    )
    return CHANGED


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
