"""Evaluate every point of a Spec and write its JSON Lines data file."""

import json
import math
import platform
import warnings
from dataclasses import dataclass
from pathlib import Path

import mpmath
import numpy as np
import scipy

from .exact import json_number, to_f64
from .spec import OVERFLOW, POLE, PRECISIONS, UNSTABLE, Spec

SCHEMA = 1
OUTPUT = "y"  # scalar functions have a single output named "y"


@dataclass(frozen=True)
class Summary:
    """Counts printed after generating one dataset."""

    name: str
    points: int
    tags: dict[str, int]


def mpmath_value(spec: Spec, args: tuple[float, ...]) -> tuple[float, bool]:
    """Evaluate at every precision in PRECISIONS (inputs are exact doubles).

    Returns the correctly rounded double and whether all precisions agree.
    """
    values = []
    for dps in PRECISIONS:
        with mpmath.workdps(dps):
            values.append(to_f64(spec.mpmath_fn(*[mpmath.mpf(a) for a in args])))
    first = values[0]
    stable = all(v == first or (math.isnan(v) and math.isnan(first)) for v in values)
    return first, stable


def scipy_value(spec: Spec, args: tuple[float, ...]) -> float:
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spec.scipy_fn(*args))


def header(spec: Spec, count: int) -> dict:
    py = platform.python_version_tuple()
    return {
        "schema": SCHEMA,
        "function": spec.name,
        "description": spec.description,
        "inputs": list(spec.inputs),
        "outputs": [OUTPUT],
        "points": count,
        "tolerance": {"rel": spec.tolerance, "reason": spec.tolerance_reason},
        "grid": spec.grid_description,
        "mpmath": {"version": mpmath.__version__, "dps": list(PRECISIONS)},
        "scipy": {"version": scipy.__version__, "function": spec.scipy_name},
        "numpy": np.__version__,
        "python": f"{py[0]}.{py[1]}",
    }


def record(spec: Spec, args: tuple[float, ...]) -> dict:
    tags = list(spec.input_tags(*args))
    if POLE in tags:
        exact = math.nan
    else:
        exact, stable = mpmath_value(spec, args)
        if not stable:
            tags.append(UNSTABLE)
        if math.isinf(exact):
            tags.append(OVERFLOW)
    if math.isnan(exact) and not tags:
        raise ValueError(f"{spec.name}{args}: mpmath gave NaN at an untagged point")
    rec = {
        "in": {name: [json_number(a)] for name, a in zip(spec.inputs, args, strict=True)},
        "mpmath": {OUTPUT: [json_number(exact)]},
        "scipy": {OUTPUT: [json_number(scipy_value(spec, args))]},
    }
    if tags:
        rec["tags"] = sorted(tags)
    return rec


def generate(spec: Spec, data_dir: Path) -> Summary:
    grid = list(spec.grid())
    records = [record(spec, args) for args in grid]
    lines = [header(spec, len(records)), *records]
    text = "".join(json.dumps(line, separators=(",", ":")) + "\n" for line in lines)
    (data_dir / f"{spec.name}.jsonl").write_text(text)
    tags: dict[str, int] = {}
    for rec in records:
        for tag in rec.get("tags", []):
            tags[tag] = tags.get(tag, 0) + 1
    return Summary(spec.name, len(records), tags)
