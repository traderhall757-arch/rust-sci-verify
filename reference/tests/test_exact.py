import math
import sys

import mpmath
import pytest

from rsv_reference.exact import json_number, to_f64


@pytest.mark.parametrize(
    "value",
    [1.0, -2.5, 0.1, 1e-310, 5e-324, -5e-324, sys.float_info.max, 1.64e-39, 2.0**-1074],
)
def test_doubles_round_trip(value):
    with mpmath.workdps(50):
        assert to_f64(mpmath.mpf(value)) == value


def test_rounds_to_nearest_including_subnormals():
    with mpmath.workdps(60):
        # Exactly halfway between 0 and the smallest subnormal: ties to even (0).
        assert to_f64(mpmath.mpf(2) ** -1075) == 0.0
        # Just above halfway rounds up.
        assert to_f64(mpmath.mpf(2) ** -1075 * (1 + mpmath.mpf(2) ** -50)) == 5e-324
        # 1/3 rounds like Python's correctly rounded division.
        assert to_f64(mpmath.mpf(1) / 3) == 1 / 3


def test_overflow_underflow_and_specials():
    with mpmath.workdps(50):
        assert to_f64(mpmath.mpf(10) ** 400) == math.inf
        assert to_f64(-(mpmath.mpf(2) ** 1024)) == -math.inf
        assert to_f64(mpmath.mpf(10) ** -400) == 0.0
        assert math.copysign(1, to_f64(-(mpmath.mpf(10) ** -400))) == -1
        assert to_f64(mpmath.mpf(0)) == 0.0
        assert math.isnan(to_f64(mpmath.nan))
        assert to_f64(mpmath.inf) == math.inf
        assert to_f64(-mpmath.inf) == -math.inf
        # The largest double plus half an ulp rounds to infinity (ties to even).
        assert to_f64(mpmath.mpf(2) ** 1024 - mpmath.mpf(2) ** 970) == math.inf


def test_json_number():
    assert json_number(math.nan) == "nan"
    assert json_number(math.inf) == "inf"
    assert json_number(-math.inf) == "-inf"
    assert json_number(0.25) == 0.25
