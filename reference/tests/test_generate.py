import math

import mpmath
import numpy as np

from rsv_reference import functions
from rsv_reference.generate import check_numpy_dispatch, generate, record
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


def test_numpy_runs_the_same_simd_code_on_every_machine():
    # The package pins NumPy below AVX-512 so SciPy values do not depend on the CPU.
    check_numpy_dispatch()
    import numpy as np

    # A point where NumPy's AVX-512 exp and the C library's differ in the last bit.
    assert float(np.exp(np.float64(-9.5))) == math.exp(-9.5)


def test_vector_inputs_and_several_outputs():
    # A sample (a tuple in the grid) is written as a list and reaches mpmath as
    # a list of mpf and SciPy as a list of floats; each output gets its own key.
    spec = _spec(
        inputs=("x", "c"),
        outputs=("total", "count"),
        grid=lambda: [((0.5, 0.25, 1e300), 2.0)],
        mpmath_fn=lambda x, c: (mpmath.fsum(x) * c, mpmath.mpf(len(x))),
        scipy_fn=lambda x, c: (sum(x) * c, float(len(x))),
    )
    rec = record(spec, ((0.5, 0.25, 1e300), 2.0))
    assert rec["in"] == {"x": [0.5, 0.25, 1e300], "c": [2.0]}
    assert rec["mpmath"] == {"total": [2e300], "count": [3.0]}
    assert rec["scipy"] == {"total": [2e300], "count": [3.0]}
    assert "tags" not in rec


def test_any_infinite_output_is_tagged_overflow():
    spec = _spec(
        outputs=("a", "b"),
        mpmath_fn=lambda x: (mpmath.exp(x), x),
        scipy_fn=lambda x: (math.inf, x),
    )
    rec = record(spec, (1000.0,))
    assert rec["tags"] == [OVERFLOW]
    assert rec["mpmath"] == {"a": ["inf"], "b": [1000.0]}


def test_wrong_number_of_outputs_is_rejected():
    spec = _spec(outputs=("a", "b"), mpmath_fn=lambda x: (x, x, x))
    try:
        record(spec, (1.0,))
    except ValueError as e:
        assert "3 outputs" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_openblas_runs_the_same_kernels_on_every_machine():
    # The package pins OpenBLAS to its Haswell kernels. At this f_oneway point
    # SciPy's statistic differs in the last bit with the SkylakeX kernels.
    check_numpy_dispatch()
    import scipy.stats as st

    r = st.f_oneway([99.5, 100.419], [101.338, 101.257])
    assert float(r.statistic) == 8.413571732372173
