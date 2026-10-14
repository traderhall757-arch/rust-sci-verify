#!/usr/bin/env python3
"""Summarise a cargo-mutants run as markdown (for the CI job summary).

Mutation score = caught / (caught + missed + timeout). Unviable mutants (ones
that do not compile) are excluded: they say nothing about the tests.

The rule (approved by the owner, see docs/guardrails.md):
- pull requests: no surviving (missed or timed-out) mutant in changed code,
  except ones excluded in .cargo/mutants.toml with a `# reason:` comment;
- weekly full run: score at least 95%.

Usage:
  scripts/mutation_score.py MUTANTS_OUT [--max-survivors N] [--min-score PCT]
  scripts/mutation_score.py --check-config .cargo/mutants.toml
"""

import argparse
import json
import re
import sys
import tomllib
from pathlib import Path


def summary(outcomes: dict) -> str:
    caught, missed, timeout = outcomes["caught"], outcomes["missed"], outcomes["timeout"]
    judged = caught + missed + timeout
    score = f"{100 * caught / judged:.1f}%" if judged else "n/a (no testable mutants)"
    lines = [
        "## Mutation testing",
        "",
        f"**Mutation score: {score}** ({caught} caught of {judged} testable mutants)",
        "",
        "| Caught | Missed | Timeout | Unviable (did not compile) |",
        "|---|---|---|---|",
        f"| {caught} | {missed} | {timeout} | {outcomes['unviable']} |",
    ]
    survivors = [
        o["scenario"]["Mutant"]["name"]
        for o in outcomes["outcomes"]
        if isinstance(o["scenario"], dict) and o["summary"] in {"MissedMutant", "Timeout"}
    ]
    if survivors:
        lines += ["", "Mutants the tests did not catch:", ""]
        lines += [f"- `{s}`" for s in survivors]
    return "\n".join(lines) + "\n"


def violations(outcomes: dict, max_survivors: int | None, min_score: float | None) -> list[str]:
    caught = outcomes["caught"]
    survivors = outcomes["missed"] + outcomes["timeout"]
    out = []
    if max_survivors is not None and survivors > max_survivors:
        out.append(f"{survivors} surviving mutant(s); at most {max_survivors} allowed")
    judged = caught + survivors
    if min_score is not None and judged and 100 * caught / judged < min_score:
        out.append(f"mutation score {100 * caught / judged:.1f}% is below {min_score:g}%")
    return out


def unexplained_exclusions(text: str) -> list[str]:
    """exclude_re entries not preceded by a `# reason:` comment line."""
    tomllib.loads(text)  # must be valid TOML
    block = re.search(r"^exclude_re\s*=\s*\[(.*?)\]\s*$", text, re.MULTILINE | re.DOTALL)
    if not block:
        return []
    out, previous = [], ""
    for line in block.group(1).splitlines():
        stripped = line.strip()
        if stripped.startswith(('"', "'")) and not previous.startswith("# reason:"):
            out.append(stripped.rstrip(","))
        if stripped:
            previous = stripped
    return out


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mutants_out", nargs="?", default="mutants.out")
    parser.add_argument("--max-survivors", type=int)
    parser.add_argument("--min-score", type=float)
    parser.add_argument("--check-config", type=Path)
    args = parser.parse_args(argv)
    if args.check_config:
        bad = unexplained_exclusions(args.check_config.read_text())
        for entry in bad:
            print(f"{args.check_config}: exclusion {entry} has no '# reason:' comment")
        print(f"mutants config check: {len(bad)} problem(s)")
        return 1 if bad else 0
    path = Path(args.mutants_out) / "outcomes.json"
    if not path.is_file():
        print(f"{path} not found: cargo-mutants did not complete", file=sys.stderr)
        return 1
    outcomes = json.loads(path.read_text())
    print(summary(outcomes), end="")
    problems = violations(outcomes, args.max_survivors, args.min_score)
    for p in problems:
        print(f"\n**FAIL:** {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
