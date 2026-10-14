import math

import mpmath
import numpy as np

from rsv_reference import functions
from rsv_reference.generate import generate, record
from rsv_reference.spec import OVERFLOW, POLE, Spec


def _spec(**kw):
    base = {
        "name": "t",
        "description": "test",
        "inputs": ("x",),
        "tolerance": 1e-15,
        "tolerance_reason": "test",
        "grid_description": "test",
        "grid": lambda: [(1.0,)],
        "mpmath_fn": mpmath.exp,
        "scipy_name": "math.exp",
        "scipy_fn": np.exp,
    }
    base.update(kw)
    return Spec(**base)


def test_overflow_is_tagged_by_rule():
    rec = record(_spec(), (1000.0,))
    assert rec["tags"] == [OVERFLOW]
    assert rec["mpmath"]["y"] == ["inf"]


def test_pole_is_tagged_from_inputs_and_not_evaluated():
    spec = _spec(
        mpmath_fn=mpmath.gamma,
        scipy_fn=lambda x: math.nan,
        input_tags=lambda x: [POLE] if x <= 0 and x == int(x) else [],
    )
    rec = record(spec, (-6.0,))
    assert rec["tags"] == [POLE]
    assert rec["mpmath"]["y"] == ["nan"]


def test_untagged_nan_is_rejected():
    spec = _spec(mpmath_fn=lambda x: mpmath.nan)
    try:
        record(spec, (1.0,))
    except ValueError as e:
        assert "untagged" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_unstable_reference_is_tagged():
    # Value depends on working precision: 1 + 10^-60 rounds differently only
    # when it is representable, i.e. at 100 digits but not at 50.
    spec = _spec(mpmath_fn=lambda x: (x + mpmath.mpf(10) ** -60 - x) * 10**60)
    rec = record(spec, (1.0,))
    assert "unstable-reference" in rec["tags"]


def test_generate_writes_header_and_points(tmp_path):
    summary = generate(_spec(grid=lambda: [(0.5,), (1000.0,)]), tmp_path)
    lines = (tmp_path / "t.jsonl").read_text().splitlines()
    assert len(lines) == 3
    assert '"schema":1' in lines[0]
    assert summary.points == 2 and summary.tags == {OVERFLOW: 1}


def test_erfcx_asymptotic_branch_matches_direct_formula():
    # At x = 1e6 both branches are valid; they must agree to working precision.
    with mpmath.workdps(50):
        x = mpmath.mpf(10) ** 6
        direct = mpmath.exp(x * x) * mpmath.erfc(x)
        series = functions._mp_erfcx(x * (1 + mpmath.mpf(10) ** -45))  # forces series branch
        assert abs(series / direct - 1) < mpmath.mpf(10) ** -40
