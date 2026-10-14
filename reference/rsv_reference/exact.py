"""Exact conversions between mpmath numbers and IEEE-754 doubles."""

import math

import mpmath

# 2**1024 is the first power of two that no finite double reaches.
_MAX_EXP2 = 1024
# Below 2**-1076 every value rounds to zero (half of the smallest subnormal).
_MIN_EXP2 = -1076


def to_f64(x: mpmath.mpf) -> float:
    """Round an mpmath real to the nearest double (ties to even).

    `float(mpf)` is not used because it may round twice for subnormal
    results. Python's integer true division is correctly rounded, so the
    exact binary value man * 2**exp is converted through it.
    """
    if mpmath.isnan(x):
        return math.nan
    if mpmath.isinf(x):
        return math.inf if x > 0 else -math.inf
    if x == 0:
        return 0.0
    # mpmath stores (sign, mantissa, exponent, bitcount); the mantissa is |x|'s.
    neg, man, exp, _ = x._mpf_
    man, exp = int(man), int(exp)
    top = exp + man.bit_length()  # |x| lies in [2**(top-1), 2**top)
    sign = -1.0 if neg else 1.0
    if top > _MAX_EXP2 + 1:
        return sign * math.inf
    if top < _MIN_EXP2:
        return sign * 0.0
    try:
        if exp >= 0:
            return sign * float(man << exp)
        return sign * (man / (1 << -exp))
    except OverflowError:
        return sign * math.inf


def json_number(v: float) -> float | str:
    """Encode a double for JSON: shortest round-trip decimal, or a string
    for the non-finite values JSON cannot represent."""
    if math.isnan(v):
        return "nan"
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return v
