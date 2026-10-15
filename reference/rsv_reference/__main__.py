"""Regenerate reference data: `python -m rsv_reference [function ...]`.

Prints one summary line per dataset; the data itself goes to reference/data/.
"""

import sys
from pathlib import Path

from .functions import ALL
from .generate import check_numpy_dispatch, generate

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def main(argv: list[str]) -> int:
    check_numpy_dispatch()
    wanted = set(argv)
    unknown = wanted - {s.name for s in ALL}
    if unknown:
        print(f"unknown function(s): {', '.join(sorted(unknown))}", file=sys.stderr)
        return 2
    DATA_DIR.mkdir(exist_ok=True)
    for spec in ALL:
        if wanted and spec.name not in wanted:
            continue
        s = generate(spec, DATA_DIR)
        tags = ", ".join(f"{n} {t}" for t, n in sorted(s.tags.items())) or "no tags"
        print(f"{s.name}: {s.points} points ({tags})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
