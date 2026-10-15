"""Regenerate reference data: `python -m rsv_reference [function ...]`.

Prints one summary line per dataset; the data itself goes to reference/data/.
(This lives outside __main__.py so the worker processes can import it.)

Every dataset is generated in its own fresh process, in parallel on all CPUs.
A dataset therefore cannot depend on state left behind by another one, which
is what lets CI regenerate only the datasets a pull request can have changed
(scripts/select_datasets.py).
"""

import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

from .functions import ALL
from .generate import Summary, generate

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _generate(name: str, data_dir: Path) -> tuple[Summary, float]:
    (spec,) = [s for s in ALL if s.name == name]
    start = time.perf_counter()
    summary = generate(spec, data_dir)
    return summary, time.perf_counter() - start


def main(argv: list[str], data_dir: Path = DATA_DIR) -> int:
    wanted = set(argv)
    unknown = wanted - {s.name for s in ALL}
    if unknown:
        print(f"unknown function(s): {', '.join(sorted(unknown))}", file=sys.stderr)
        return 2
    names = [s.name for s in ALL if not wanted or s.name in wanted]
    data_dir.mkdir(exist_ok=True)
    with ProcessPoolExecutor(
        max_workers=min(len(names), os.cpu_count() or 1),
        mp_context=multiprocessing.get_context("spawn"),
        max_tasks_per_child=1,  # a fresh process per dataset
    ) as pool:
        for s, seconds in pool.map(_generate, names, repeat(data_dir)):
            tags = ", ".join(f"{n} {t}" for t, n in sorted(s.tags.items())) or "no tags"
            print(f"{s.name}: {s.points} points ({tags}) in {seconds:.1f} s", flush=True)
    return 0
