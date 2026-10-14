#!/usr/bin/env python3
"""Summarise a cargo-mutants run as markdown (for the CI job summary).

Mutation score = caught / (caught + missed + timeout). Unviable mutants (ones
that do not compile) are excluded: they say nothing about the tests.
Report-only for now; the owner sets a minimum after Milestone 1b.

Usage: scripts/mutation_score.py [mutants.out]
"""

import json
import sys
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


def main(argv: list[str]) -> int:
    path = Path(argv[0] if argv else "mutants.out") / "outcomes.json"
    if not path.is_file():
        print(f"{path} not found: cargo-mutants did not complete", file=sys.stderr)
        return 1
    print(summary(json.loads(path.read_text())), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
