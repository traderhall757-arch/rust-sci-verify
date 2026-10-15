"""Evaluate every point of a Spec and write its JSON Lines data file."""

import json
import math
import os
import platform
import warnings
from dataclasses import dataclass
from pathlib import Path

import mpmath
import numpy as np
import scipy

from . import NUMPY_DISABLED_TARGETS, OPENBLAS_CORETYPE
from .exact import json_number, to_f64
from .spec import OVERFLOW, POLE, PRECISIONS, UNSTABLE, Spec


def check_numpy_dispatch() -> None:
    """Fail if NumPy runs AVX-512 code or OpenBLAS is not pinned (see the
    package's notes).

    NumPy reads NPY_DISABLE_CPU_FEATURES only when it is first imported, so
    an earlier import (for example by a test runner) defeats the pin.
    """
    from numpy._core._multiarray_umath import __cpu_features__

    active = [t for t in NUMPY_DISABLED_TARGETS.split() if __cpu_features__.get(t)]
    if active:
        raise RuntimeError(
            f"NumPy AVX-512 targets active ({', '.join(active)}): reference values would "
            "depend on the CPU. Import rsv_reference before numpy."
        )
    if os.environ.get("OPENBLAS_CORETYPE") != OPENBLAS_CORETYPE:
        raise RuntimeError(
            f"OPENBLAS_CORETYPE is not {OPENBLAS_CORETYPE}: reference values would depend on "
            "the CPU. Import rsv_reference before numpy and scipy."
        )


SCHEMA = 1


@dataclass(frozen=True)
class Summary:
    """Counts printed after generating one dataset."""

    name: str
    points: int
    tags: dict[str, int]


Arg = float | tuple[float, ...]


def _mp_arg(a: Arg) -> mpmath.mpf | list[mpmath.mpf]:
    return [mpmath.mpf(v) for v in a] if isinstance(a, tuple) else mpmath.mpf(a)


def _outputs(spec: Spec, value) -> tuple:
    """The function's value as a tuple with one entry per output."""
    if len(spec.outputs) == 1:
        return (value,)
    value = tuple(value)
    if len(value) != len(spec.outputs):
        raise ValueError(f"{spec.name}: {len(value)} outputs, expected {len(spec.outputs)}")
    return value


def _same(a: float, b: float) -> bool:
    return a == b or (math.isnan(a) and math.isnan(b))


def mpmath_value(spec: Spec, args: tuple[Arg, ...]) -> tuple[tuple[float, ...], bool]:
    """Evaluate at every precision in PRECISIONS (inputs are exact doubles).

    Returns the correctly rounded doubles (one per output) and whether all
    precisions agree on every output.
    """
    values = []
    for dps in PRECISIONS:
        with mpmath.workdps(dps):
            out = _outputs(spec, spec.mpmath_fn(*[_mp_arg(a) for a in args]))
            values.append(tuple(to_f64(v) for v in out))
    first = values[0]
    stable = all(_same(a, b) for v in values for a, b in zip(v, first, strict=True))
    return first, stable


def scipy_value(spec: Spec, args: tuple[Arg, ...]) -> tuple[float, ...]:
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sci_args = [list(a) if isinstance(a, tuple) else a for a in args]
        return tuple(float(v) for v in _outputs(spec, spec.scipy_fn(*sci_args)))


def header(spec: Spec, count: int) -> dict:
    py = platform.python_version_tuple()
    return {
        "schema": SCHEMA,
        "function": spec.name,
        "description": spec.description,
        "inputs": list(spec.inputs),
        "outputs": list(spec.outputs),
        "points": count,
        "tolerance": {"rel": spec.tolerance, "reason": spec.tolerance_reason},
        "grid": spec.grid_description,
        "mpmath": {"version": mpmath.__version__, "dps": list(PRECISIONS)},
        "scipy": {"version": scipy.__version__, "function": spec.scipy_name},
        "numpy": np.__version__,
        "python": f"{py[0]}.{py[1]}",
    }


def _json_list(a: Arg) -> list:
    return [json_number(v) for v in a] if isinstance(a, tuple) else [json_number(a)]


def record(spec: Spec, args: tuple[Arg, ...]) -> dict:
    tags = list(spec.input_tags(*args))
    if POLE in tags:
        exact = tuple(math.nan for _ in spec.outputs)
    else:
        exact, stable = mpmath_value(spec, args)
        if not stable:
            tags.append(UNSTABLE)
        if any(math.isinf(v) for v in exact):
            tags.append(OVERFLOW)
    if any(math.isnan(v) for v in exact) and not tags:
        raise ValueError(f"{spec.name}{args}: mpmath gave NaN at an untagged point")
    scipy_out = scipy_value(spec, args)
    rec = {
        "in": {name: _json_list(a) for name, a in zip(spec.inputs, args, strict=True)},
        "mpmath": {o: [json_number(v)] for o, v in zip(spec.outputs, exact, strict=True)},
        "scipy": {o: [json_number(v)] for o, v in zip(spec.outputs, scipy_out, strict=True)},
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
