"""The registry of reference datasets. One entry per mathematical function."""

import itertools
import math
from fractions import Fraction

import mpmath
import scipy.special as sc
import scipy.stats as st

from .spec import POLE, Spec


def _decades() -> list[float]:
    # 1e-300 .. 1e300 every 10 decades. float("1eK") is correctly rounded, so
    # the grid does not depend on the platform's pow().
    return [float(f"1e{k}") for k in range(-300, 301, 10)]


def _erfcx_grid() -> list[tuple[float]]:
    # Dense uniform part, step 0.25 (exact in binary), covering the
    # transition from the series/continued-fraction regions.
    dense = [i / 4 for i in range(-120, 121)]
    decades = _decades()
    points = sorted(set(dense + decades + [-d for d in decades]))
    return [(x,) for x in points]


def _mp_erfcx(x: mpmath.mpf) -> mpmath.mpf:
    if x < 10**6:
        return mpmath.exp(x * x) * mpmath.erfc(x)
    # mpmath's erfc overflows internally for x > ~1e154, so large x uses the
    # asymptotic series erfcx(x) = 1/(x sqrt(pi)) * sum (-1)^n (2n-1)!! / (2x^2)^n
    # (Abramowitz & Stegun 7.1.23). For real x > 0 the series is enveloping:
    # the error is smaller than the first omitted term, which we drive below
    # the working precision.
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    total, term, n = mpmath.mpf(1), mpmath.mpf(1), 0
    while abs(term) > eps:
        n += 1
        term *= -(2 * n - 1) / (2 * x * x)
        total += term
    return total / (x * mpmath.sqrt(mpmath.pi))


ERFCX = Spec(
    name="erfcx",
    description="Scaled complementary error function erfcx(x) = exp(x^2) erfc(x)",
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "For negative x the value is about 2*exp(x^2), and merely rounding x^2 to a "
        "double costs up to x^2 * 2^-53 relative error (about 8e-14 at x = -26.6, "
        "the overflow edge); 1e-13 accepts that unavoidable loss and nothing more"
    ),
    grid_description=(
        "x = -30..30 step 0.25, plus +-10^k for k = -300..300 step 10; "
        "negative x below about -26.6 overflows"
    ),
    grid=_erfcx_grid,
    mpmath_fn=_mp_erfcx,
    scipy_name="scipy.special.erfcx",
    scipy_fn=sc.erfcx,
)

# --- gamma -----------------------------------------------------------------


def _gamma_grid() -> list[tuple[float]]:
    # x = -30..180 step 0.25 (exact in binary): every pole from 0 to -30, the
    # half-integers between them, and the overflow edge near 171.62; plus
    # +-10^-k for tiny |x|, where gamma(x) is about 1/x.
    dense = [i / 4 for i in range(-120, 721)]
    tiny = [d for d in _decades() if d < 1e-5]
    points = sorted(set(dense + tiny + [-d for d in tiny]))
    return [(x,) for x in points]


def _gamma_tags(x: float) -> list[str]:
    return [POLE] if x <= 0 and x == int(x) else []


GAMMA = Spec(
    name="gamma",
    description="Gamma function gamma(x) for real x",
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "gamma(x) grows like ((x + g - 1/2)/e)^(x - 1/2); rounding the base to a double is "
        "amplified by the exponent, up to 171 * 2^-53 (about 2e-14) at the overflow edge "
        "x = 171.6, and by the same amount through the reflection formula for negative x; "
        "1e-13 allows that plus a few ulps of further arithmetic"
    ),
    grid_description=(
        "x = -30..180 step 0.25 (poles at 0, -1, ..., -30; overflow above about 171.62), "
        "plus +-10^k for k = -300..-10 step 10"
    ),
    grid=_gamma_grid,
    mpmath_fn=mpmath.gamma,
    scipy_name="scipy.special.gamma",
    scipy_fn=sc.gamma,
    input_tags=_gamma_tags,
)


# --- erfc ------------------------------------------------------------------


def _erfc_grid() -> list[tuple[float]]:
    # Dense part through the whole range where erfc is not 0 or 2 in double
    # precision (it underflows completely near x = 27.3), plus decades.
    dense = [i / 4 for i in range(-40, 121)]
    decades = _decades()
    points = sorted(set(dense + decades + [-d for d in decades]))
    return [(x,) for x in points]


def _mp_erfc(x: mpmath.mpf) -> mpmath.mpf:
    if x < 10**6:
        return mpmath.erfc(x)
    # mpmath's erfc overflows internally for x > ~1e154; use the tested
    # asymptotic series of erfcx instead (mpmath numbers do not underflow).
    return _mp_erfcx(x) * mpmath.exp(-x * x)


ERFC = Spec(
    name="erfc",
    description="Complementary error function erfc(x) = 1 - erf(x)",
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "For positive x the value is about exp(-x^2)/(x sqrt(pi)), and merely rounding x^2 "
        "to a double costs up to x^2 * 2^-53 relative error (about 8e-14 at x = 26.5, where "
        "the result becomes subnormal); 1e-13 accepts that unavoidable loss and nothing more"
    ),
    grid_description=(
        "x = -10..30 step 0.25, plus +-10^k for k = -300..300 step 10; "
        "results are subnormal from about x = 26.6 and zero beyond about 27.3"
    ),
    grid=_erfc_grid,
    mpmath_fn=_mp_erfc,
    scipy_name="scipy.special.erfc",
    scipy_fn=sc.erfc,
)


# --- gammainc --------------------------------------------------------------

_GAMMAINC_A = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0)
_GAMMAINC_X = (
    0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
    100.0, 200.0, 500.0, 1000.0, 2000.0,
)  # fmt: skip


def _gammainc_grid() -> list[tuple[float, float]]:
    # Every a is also an x, so the grid includes x = a, where P is close to
    # 1/2 and the series / continued-fraction transition lies.
    return sorted(itertools.product(_GAMMAINC_A, _GAMMAINC_X))


def _mp_gammainc(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    return mpmath.gammainc(a, 0, x, regularized=True)


GAMMAINC = Spec(
    name="gammainc",
    description="Regularized lower incomplete gamma function P(a, x) = gamma(a, x) / Gamma(a)",
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "P is the prefix x^a e^-x / Gamma(a + 1) times a well-conditioned series or continued "
        "fraction. On this grid the logarithms of the prefix's factors reach 1.6e4 in "
        "magnitude while their sum stays above -745; forming them separately in double "
        "precision loses up to 1.6e4 * 2^-53 (about 2e-12), and combining them before "
        "rounding (Temme; Boost) loses far less. 1e-12 demands the careful approach and is "
        "3x the worst SciPy error on this grid (3.5e-13 at a = 1000, x = 500)"
    ),
    grid_description=(
        "a in {0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000} x "
        "x in {0, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000}"
    ),
    grid=_gammainc_grid,
    mpmath_fn=_mp_gammainc,
    scipy_name="scipy.special.gammainc",
    scipy_fn=sc.gammainc,
)


# --- betainc ---------------------------------------------------------------

_BETAINC_AB = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000.0, 10000.0)
_BETAINC_X = (0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)


def _betainc_grid() -> list[tuple[float, float, float]]:
    # Every (a, b, x) triple, plus x = a/(a+b) (the mode region, where the
    # result is far from 0 and 1 even for large a and b).
    triples = set(itertools.product(_BETAINC_AB, _BETAINC_AB, _BETAINC_X))
    triples |= {(a, b, a / (a + b)) for a in _BETAINC_AB for b in _BETAINC_AB}
    return sorted(triples)


def _mp_betainc(a: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """I_x(a, b) for 0 < x < 1 by the hypergeometric series DLMF 8.17.8:

        I_x(a, b) = x^a (1-x)^b / (a B(a, b)) * sum_n (a+b)_n / (a+1)_n x^n,

    used for x <= (a+1)/(a+b+2), where the term ratio (a+b+n) x / (a+1+n) is
    below 1 from the start, and the symmetry I_x(a, b) = 1 - I_(1-x)(b, a)
    otherwise. All terms are positive, so the sum has no cancellation.

    mpmath.betainc is not used because it takes minutes per point near the
    mode for a, b around 1e4. Tested against mpmath.betainc and an exact
    binomial sum in tests/test_functions.py.
    """
    y = mpmath.fsub(1, x, exact=True)
    if x > (a + 1) / (a + b + 2):
        return 1 - _mp_betainc(b, a, y)
    eps = mpmath.mpf(2) ** (-mpmath.mp.prec - 10)
    total = term = mpmath.mpf(1)
    n = 0
    while True:
        ratio = (a + b + n) * x / (a + 1 + n)
        term *= ratio
        total += term
        n += 1
        # The ratio decreases with n when b >= 1 and stays below x when
        # b < 1, which bounds the remaining tail geometrically.
        bound = ratio if b >= 1 else x
        if term * bound / (1 - bound) < eps * total:
            break
    return x**a * y**b / (a * mpmath.beta(a, b)) * total


BETAINC = Spec(
    name="betainc",
    description="Regularized incomplete beta function I_x(a, b) = B(x; a, b) / B(a, b)",
    inputs=("a", "b", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "I is the prefix x^a (1-x)^b / (a B(a, b)) times a well-conditioned series or "
        "continued fraction. On this grid the logarithms of the prefix's factors reach "
        "2.8e4 in magnitude (a = b = 1e4) while their sum stays above -745; forming them "
        "separately in double precision loses up to 2.8e4 * 2^-53 (about 3e-12), and "
        "combining them before rounding (Boost's Lanczos-based prefix) loses far less. "
        "1e-12 demands the careful approach and is 6x the worst SciPy error on this grid "
        "(1.6e-13)"
    ),
    grid_description=(
        "a, b in {0.1, 0.5, 1, 2, 5, 10, 100, 1000, 10000} x "
        "x in {0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99}, plus x = a/(a+b)"
    ),
    grid=_betainc_grid,
    mpmath_fn=_mp_betainc,
    scipy_name="scipy.special.betainc",
    scipy_fn=sc.betainc,
)


# --- betainc_large -----------------------------------------------------------

_BETAINC_LARGE_AB = (3001.0, 5000.0, 20000.0, 100000.0)
_BETAINC_LARGE_K = (-10, -5, -2, -1, 0, 1, 2, 5, 10)


def _betainc_large_grid() -> list[tuple[float, float, float]]:
    # Both parameters above 3000, where puruspe switches to Gauss-Legendre
    # quadrature. x = mean + k standard deviations of the beta distribution,
    # so every row covers both tails (I near 0 and near 1) and the centre.
    points = set()
    for a, b in itertools.product(_BETAINC_LARGE_AB, repeat=2):
        mean = a / (a + b)
        sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
        points |= {(a, b, mean + k * sd) for k in _BETAINC_LARGE_K}
    return sorted(points)


BETAINC_LARGE = Spec(
    name="betainc_large",
    description=(
        "Regularized incomplete beta function I_x(a, b) for large a and b "
        "(same function as betainc)"
    ),
    inputs=("a", "b", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "Same reasoning as betainc: the logarithms of the prefix's factors reach "
        "2.8e5 in magnitude here (a = b = 1e5), so forming them separately in double "
        "precision loses up to 2.8e5 * 2^-53 (about 3e-11), while combining them before "
        "rounding does not. 1e-12 is 7x the worst SciPy error on this grid (1.5e-13)"
    ),
    grid_description=(
        "a, b in {3001, 5000, 2e4, 1e5} x x = a/(a+b) + k sd, "
        "k in {-10, -5, -2, -1, 0, 1, 2, 5, 10}, sd the standard deviation of Beta(a, b)"
    ),
    grid=_betainc_large_grid,
    mpmath_fn=_mp_betainc,
    scipy_name="scipy.special.betainc",
    scipy_fn=sc.betainc,
)


# --- betainc_switch ----------------------------------------------------------

_BETAINC_SWITCH_A = (1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
_BETAINC_SWITCH_B = (0.01, 0.02, 0.05, 0.2, 0.5)
_BETAINC_SWITCH_Y = (1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
_BETAINC_SWITCH_AB = (10.0, 100.0, 1000.0, 2900.0)
# Worst points of the 2026-10-14 random spot-check (finding 001), as quoted there.
_BETAINC_SWITCH_SPOT = (
    (73256.2, 0.01866, 0.9999945374),
    (58386.6, 0.01965, 0.9999866779),
    (97698.9, 0.0166, 0.9999973151),
    (7457.77, 0.02365, 0.9998657733),
    (15.766, 0.01229, 0.9565398173),
)


def _betainc_switch_grid() -> list[tuple[float, float, float]]:
    # Where the choice between the direct form I_x(a, b) and the complement
    # 1 - I_(1-x)(b, a) decides the accuracy (finding 001):
    # - small b with x close to 1, where I is still far from 1: x = 1 - y;
    # - b > 2, halfway between puruspe 0.4.4's switch (a+1)/(a+2b) and
    #   Numerical Recipes' (a+1)/(a+b+2).
    # Left out: y < 1e-2 with a * y >= 2 (only a = 1e5, y = 1e-4), where x is
    # below every candidate switch point, so all use the direct fraction,
    # and the reference series needs about 1/y terms (minutes per point).
    small_b = {
        (a, b, 1 - y)
        for a, b, y in itertools.product(_BETAINC_SWITCH_A, _BETAINC_SWITCH_B, _BETAINC_SWITCH_Y)
        if y >= 1e-2 or a * y < 2
    }
    between = {
        (a, b, ((a + 1) / (a + 2 * b) + (a + 1) / (a + b + 2)) / 2)
        for a, b in itertools.product(_BETAINC_SWITCH_AB, repeat=2)
    }
    # Points already in betainc are left out.
    return sorted((small_b | between | set(_BETAINC_SWITCH_SPOT)) - set(_betainc_grid()))


BETAINC_SWITCH = Spec(
    name="betainc_switch",
    description=(
        "Regularized incomplete beta function I_x(a, b) where the choice between "
        "I_x(a, b) and 1 - I_(1-x)(b, a) matters (same function as betainc)"
    ),
    inputs=("a", "b", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "Same reasoning as betainc: the logarithms of the prefix's factors reach 1.1e6 in "
        "magnitude here (ln Gamma(1e5)), so forming them separately in double precision loses "
        "up to about 1e-10, and computing a small I as 1 - I_(1-x)(b, a) magnifies that by "
        "(1 - I)/I. Combining them before rounding, and computing the smaller of I and 1 - I "
        "directly, avoids both. 1e-12 is 17x the worst SciPy error on this grid (5.9e-14)"
    ),
    grid_description=(
        "a in {1, 10, 100, 1e3, 1e4, 1e5} x b in {0.01, 0.02, 0.05, 0.2, 0.5} x "
        "x = 1 - y, y in {1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10}; a, b in {10, 100, 1000, "
        "2900} with x halfway between (a+1)/(a+2b) and (a+1)/(a+b+2); 5 points from a "
        "random spot-check (finding 001); without a = 1e5 with y = 1e-4 and points in betainc"
    ),
    grid=_betainc_switch_grid,
    mpmath_fn=_mp_betainc,
    scipy_name="scipy.special.betainc",
    scipy_fn=sc.betainc,
)


# --- kv --------------------------------------------------------------------

_KV_NU = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0)
_KV_X = (
    0.01, 0.0631, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 700.0,
)  # fmt: skip


def _kv_grid() -> list[tuple[float, float]]:
    # Integer orders, so the same data checks integer-order implementations.
    # x = 0.0631 is the pilot point where SciPy returns inf for nu = 100.
    return sorted(itertools.product(_KV_NU, _KV_X))


KV = Spec(
    name="kv",
    description="Modified Bessel function of the second kind K_nu(x), integer order",
    inputs=("nu", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "For small x, K_nu(x) is about Gamma(nu)/2 * (2/x)^nu; rounding the base to a double "
        "is amplified by the exponent, up to 100 * 2^-53 (about 1e-14) at nu = 100; for large "
        "x it is about sqrt(pi/(2x)) e^-x with x exact. 1e-13 allows that plus a few ulps of "
        "further arithmetic"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100} x x in {0.01, 0.0631, 0.1, 0.25, 0.5, 1, 2, 5, "
        "10, 20, 50, 100, 200, 500, 700}; K_100(0.01) overflows"
    ),
    grid=_kv_grid,
    mpmath_fn=mpmath.besselk,
    scipy_name="scipy.special.kv",
    scipy_fn=sc.kv,
)


# --- yv --------------------------------------------------------------------

# Integer orders, so the same data checks integer-order implementations.
_YV_NU = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0)
_YV_X = (1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
_YV_LARGE_X = (50.0, 100.0, 200.0, 500.0, 1e3, 1e4, 1e5, 1e6)


def _yv_grid() -> list[tuple[float, float]]:
    return sorted(itertools.product(_YV_NU, _YV_X))


def _yv_large_grid() -> list[tuple[float, float]]:
    return sorted(itertools.product(_YV_NU, _YV_LARGE_X))


def yv_condition(nu: float, x: float) -> mpmath.mpf:
    """|x Y_nu'(x) / Y_nu(x)|: the relative change of Y per relative change of x.

    An error of e (relative) in any quantity proportional to x, such as the
    phase x - (2 nu + 1) pi / 4 of the large-x expansion rounded to a double,
    changes Y by about this times e.
    """
    return abs(x * mpmath.bessely(nu, x, derivative=1) / mpmath.bessely(nu, x))


YV = Spec(
    name="yv",
    description="Bessel function of the second kind Y_nu(x), integer order, x <= 20",
    inputs=("nu", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "Rounding a quantity proportional to x to a double changes Y by the condition number "
        "|x Y'/Y| times 2^-53. On this grid the condition number is at most 430 (about nu for "
        "small x; 426 at nu = 2, x = 10, near a zero), so up to 4.8e-14 is unavoidable; 1e-13 allows that plus a few ulps of further arithmetic, as for kv"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {1e-10, 1e-5, 1e-3, 0.01, 0.1, "
        "0.25, 0.5, 1, 2, 5, 10, 20}; Y overflows to -inf for large nu at small x"
    ),
    grid=_yv_grid,
    mpmath_fn=mpmath.bessely,
    scipy_name="scipy.special.yv",
    scipy_fn=sc.yv,
)

YV_LARGE = Spec(
    name="yv_large",
    description="Bessel function of the second kind Y_nu(x), integer order, x from 50 to 1e6",
    inputs=("nu", "x"),
    tolerance=5e-10,
    tolerance_reason=(
        "For large x, Y oscillates with phase about x - (2 nu + 1) pi / 4, which cannot be "
        "held exactly in a double: rounding it costs the condition number |x Y'/Y| times "
        "2^-53. The condition number grows like x and reaches 2.2e6 on this grid (nu = 1, "
        "x = 1e6), so up to 2.4e-10 is unavoidable for a method that rounds the phase; 5e-10 "
        "allows 2x that and is 19x the worst SciPy error on this grid (2.6e-11 at nu = 500, "
        "x = 1e6). It is loose for the smaller x here: the report shows each error"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {50, 100, 200, 500, 1e3, 1e4, "
        "1e5, 1e6}"
    ),
    grid=_yv_large_grid,
    mpmath_fn=mpmath.bessely,
    scipy_name="scipy.special.yv",
    scipy_fn=sc.yv,
)


# --- jv --------------------------------------------------------------------

# The same grids as yv (integer orders), so J and Y are checked at the same points.


def jv_condition(nu: float, x: float) -> mpmath.mpf:
    """|x J_nu'(x) / J_nu(x)|: the relative change of J per relative change of x."""
    return abs(x * mpmath.besselj(nu, x, derivative=1) / mpmath.besselj(nu, x))


JV = Spec(
    name="jv",
    description="Bessel function of the first kind J_nu(x), integer order, x <= 20",
    inputs=("nu", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "Rounding a quantity proportional to x to a double changes J by the condition number "
        "|x J'/J| times 2^-53. On this grid the condition number is at most 200 (about nu for "
        "small x), so up to 2.2e-14 is unavoidable; 1e-13 allows that plus a few ulps of "
        "further arithmetic, as for yv. Results below 2.2e-308 are judged relative to that "
        "number (conventions), so the underflowed J_500(1e-10) etc. need no special rule"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {1e-10, 1e-5, 1e-3, 0.01, 0.1, "
        "0.25, 0.5, 1, 2, 5, 10, 20} (the yv grid); J underflows to 0 for large nu at small x"
    ),
    grid=_yv_grid,
    mpmath_fn=mpmath.besselj,
    scipy_name="scipy.special.jv",
    scipy_fn=sc.jv,
)

JV_LARGE = Spec(
    name="jv_large",
    description="Bessel function of the first kind J_nu(x), integer order, x from 50 to 1e6",
    inputs=("nu", "x"),
    tolerance=5e-10,
    tolerance_reason=(
        "For large x, J oscillates with phase about x - (2 nu + 1) pi / 4, which cannot be "
        "held exactly in a double: rounding it costs the condition number |x J'/J| times "
        "2^-53. The condition number grows like x and reaches 2.2e6 on this grid (nu = 0, "
        "x = 1e6), so up to 2.4e-10 is unavoidable for a method that rounds the phase; 5e-10 "
        "allows 2x that, the same as yv_large; SciPy's worst on this grid is 2.0e-10 (nu = "
        "200, x = 1e4). It is loose for the smaller x here: the report shows each error"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {50, 100, 200, 500, 1e3, 1e4, "
        "1e5, 1e6} (the yv_large grid)"
    ),
    grid=_yv_large_grid,
    mpmath_fn=mpmath.besselj,
    scipy_name="scipy.special.jv",
    scipy_fn=sc.jv,
)

# --- iv --------------------------------------------------------------------

# Integer orders on the yv grids (so In can use the data), but larger x stops
# before the overflow edge (I_0(x) overflows near x = 713.99) and a little past it.
_IV_LARGE_X = (50.0, 100.0, 200.0, 500.0, 700.0, 710.0, 720.0, 1000.0)
# Non-integer orders, for implementations of real order only.
_IV_NONINT_NU = (0.25, 0.5, 1.5, 2.7, 10.5, 25.25, 100.5, 200.75)
_IV_NONINT_X = (1e-5, 0.01, 0.1, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 500.0, 700.0)


def _iv_large_grid() -> list[tuple[float, float]]:
    return sorted(itertools.product(_YV_NU, _IV_LARGE_X))


def _iv_nonint_grid() -> list[tuple[float, float]]:
    return sorted(itertools.product(_IV_NONINT_NU, _IV_NONINT_X))


def iv_condition(nu: float, x: float) -> mpmath.mpf:
    """|x I_nu'(x) / I_nu(x)|: the relative change of I per relative change of x.

    About nu for small x and about x for large x (I grows like e^x).
    """
    return abs(x * mpmath.besseli(nu, x, derivative=1) / mpmath.besseli(nu, x))


IV = Spec(
    name="iv",
    description="Modified Bessel function of the first kind I_nu(x), integer order, x <= 20",
    inputs=("nu", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "Rounding a quantity proportional to x to a double changes I by the condition number "
        "|x I'/I| times 2^-53. On this grid the condition number is at most 201 (about nu for "
        "small x), so up to 2.2e-14 is unavoidable; 1e-13 allows that plus a few ulps of "
        "further arithmetic, as for jv; SciPy's worst on this grid is 6.2e-14 (nu = 100, "
        "x = 0.1). Results below 2.2e-308 are judged relative to that number (conventions), "
        "so the underflowed I_500(1e-10) etc. need no special rule"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {1e-10, 1e-5, 1e-3, 0.01, 0.1, "
        "0.25, 0.5, 1, 2, 5, 10, 20} (the yv grid); I underflows to 0 for large nu at small x"
    ),
    grid=_yv_grid,
    mpmath_fn=mpmath.besseli,
    scipy_name="scipy.special.iv",
    scipy_fn=sc.iv,
)

IV_LARGE = Spec(
    name="iv_large",
    description="Modified Bessel function of the first kind I_nu(x), integer order, x from 50 to 1000",
    inputs=("nu", "x"),
    tolerance=2e-13,
    tolerance_reason=(
        "For large x, I grows like e^(nu eta) with nu eta about x (uniform expansion, DLMF "
        "10.41.3); rounding that exponent to a double costs the condition number |x I'/I| "
        "times 2^-53. On this grid the condition number is at most 876 at a finite point "
        "(nu = 500, x = 720), so up to 9.7e-14 is unavoidable for a method that rounds the "
        "exponent; 2e-13 allows 2x that. SciPy's worst finite result here is 1.7e-13 (nu = "
        "500, x = 710); it returns inf for I_100(720) = 7.1e307"
    ),
    grid_description=(
        "nu in {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} x x in {50, 100, 200, 500, 700, 710, "
        "720, 1000}; I overflows from x = 720 for nu <= 50 and at x = 1000 for every nu"
    ),
    grid=_iv_large_grid,
    mpmath_fn=mpmath.besseli,
    scipy_name="scipy.special.iv",
    scipy_fn=sc.iv,
)

IV_NONINT = Spec(
    name="iv_nonint",
    description="Modified Bessel function of the first kind I_nu(x), non-integer order",
    inputs=("nu", "x"),
    tolerance=2e-13,
    tolerance_reason=(
        "Same reasoning as iv and iv_large: the condition number |x I'/I| is at most 728 on "
        "this grid (nu = 200.75, x = 700), so up to 8.1e-14 is unavoidable; 2e-13 allows 2x "
        "that, the same as iv_large. SciPy's worst on this grid is 8.0e-14 (nu = 100.5, "
        "x = 0.1)"
    ),
    grid_description=(
        "nu in {0.25, 0.5, 1.5, 2.7, 10.5, 25.25, 100.5, 200.75} x x in {1e-5, 0.01, 0.1, 1, "
        "2, 5, 10, 20, 50, 100, 500, 700}; I underflows to 0 for large nu at small x"
    ),
    grid=_iv_nonint_grid,
    mpmath_fn=mpmath.besseli,
    scipy_name="scipy.special.iv",
    scipy_fn=sc.iv,
)


# --- standard normal distribution ----------------------------------------------


def _norm_x_grid() -> list[tuple[float]]:
    # x = -40..40 step 0.25 (exact in binary) covers both tails down to complete
    # underflow (the lower tail is subnormal from about -37.52 and 0 below about
    # -38.49) plus +-10^k decades for tiny and huge |x|.
    dense = [i / 4 for i in range(-160, 161)]
    decades = _decades()
    points = sorted(set(dense + decades + [-d for d in decades]))
    return [(x,) for x in points]


def _mp_norm_cdf(x: mpmath.mpf) -> mpmath.mpf:
    # Phi(x) = erfc(-x / sqrt(2)) / 2 (DLMF 7.18 / 26.2.29 of A&S), with the
    # large-argument branch of _mp_erfc for |x| >= sqrt(2) * 1e6.
    return _mp_erfc(-x / mpmath.sqrt(2)) / 2


def _mp_norm_sf(x: mpmath.mpf) -> mpmath.mpf:
    return _mp_norm_cdf(-x)


def norm_cdf_condition(x: float) -> mpmath.mpf:
    """|x phi(x) / Phi(x)|: the relative change of Phi per relative change of x.

    About x^2 in the lower tail (Phi(x) is about phi(x) / |x| there); the
    upper-tail survival function has the same condition number at -x.
    """
    x = mpmath.mpf(x)
    return abs(x * mpmath.npdf(x) / _mp_norm_cdf(x))


NORM_CDF = Spec(
    name="norm_cdf",
    description="Standard normal cumulative distribution function Phi(x) (mean 0, sd 1)",
    inputs=("x",),
    tolerance=4e-13,
    tolerance_reason=(
        "Rounding a quantity proportional to x (such as x / sqrt(2)) to a double changes Phi "
        "by the condition number |x phi(x) / Phi(x)| times 2^-53. In the lower tail the "
        "condition number is about x^2 and reaches 1408 at the last normal result "
        "(x = -37.5), so up to 1.6e-13 is unavoidable for a method that rounds such a "
        "quantity; 4e-13 allows 2x that, as for iv_large. SciPy's worst normal result here "
        "is 1.1e-13 (x = -32.25); it flushes subnormal results at x = -37.75 and -38 to 0. "
        "Results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        "x = -40..40 step 0.25, plus +-10^k for k = -300..300 step 10; Phi is subnormal "
        "from about x = -37.52, 0 below about -38.49 and rounds to 1 above about 8.3"
    ),
    grid=_norm_x_grid,
    mpmath_fn=_mp_norm_cdf,
    scipy_name="scipy.stats.norm.cdf",
    scipy_fn=st.norm.cdf,
)

NORM_SF = Spec(
    name="norm_sf",
    description="Standard normal survival function Q(x) = 1 - Phi(x) = Phi(-x) (mean 0, sd 1)",
    inputs=("x",),
    tolerance=4e-13,
    tolerance_reason=(
        "Q(x) = Phi(-x), so the reasoning of norm_cdf applies mirrored: the condition number "
        "|x phi(x) / Q(x)| is about x^2 in the upper tail and reaches 1408 at x = 37.5, so up "
        "to 1.6e-13 is unavoidable; 4e-13 allows 2x that. SciPy's worst normal result here is "
        "1.1e-13 (x = 32.25)"
    ),
    grid_description=(
        "x = -40..40 step 0.25, plus +-10^k for k = -300..300 step 10 (the norm_cdf grid); "
        "Q is subnormal from about x = 37.52, 0 above about 38.49 and rounds to 1 below "
        "about -8.3"
    ),
    grid=_norm_x_grid,
    mpmath_fn=_mp_norm_sf,
    scipy_name="scipy.stats.norm.sf",
    scipy_fn=st.norm.sf,
)


def _norm_ppf_grid() -> list[tuple[float]]:
    # Both tails (p = 10^-k and the doubles nearest 1 - 10^-k), the smallest
    # subnormal, the bulk in steps of 1/64 (exact), and points just off 1/2.
    tail = [float(f"1e-{k}") for k in range(1, 21)] + [float(f"1e-{k}") for k in range(30, 301, 10)]
    upper = [1 - float(f"1e-{k}") for k in range(1, 17)]
    bulk = [i / 64 for i in range(65)]
    near_half = [0.5 + s * float(f"1e-{k}") for k in (5, 10, 15) for s in (-1, 1)]
    extremes = [5e-324, 1e-310, 2.0**-1022]
    return [(p,) for p in sorted(set(tail + upper + bulk + near_half + extremes))]


def _norm_ppf_tags(p: float) -> list[str]:
    # The quantile is -inf at p = 0 and +inf at p = 1.
    return [POLE] if p in (0.0, 1.0) else []


def _mp_norm_ppf(p: mpmath.mpf) -> mpmath.mpf:
    """The x with Phi(x) = p, by Newton's method on log Phi (0 < p < 1).

    p is an exact double, so 1 - p is exact at any working precision above 53
    bits and the upper half uses the symmetry ppf(p) = -ppf(1 - p).
    """
    if p == mpmath.mpf(1) / 2:
        return mpmath.mpf(0)
    if p > mpmath.mpf(1) / 2:
        return -_mp_norm_ppf(1 - p)
    # Start from A&S 26.2.23 (error below 4.5e-4) and iterate
    # x <- x - (log Phi(x) - log p) Phi(x) / phi(x), which converges
    # quadratically from this start.
    t = mpmath.sqrt(-2 * mpmath.log(p))
    x = -(
        t
        - (2.515517 + 0.802853 * t + 0.010328 * t**2)
        / (1 + 1.432788 * t + 0.189269 * t**2 + 0.001308 * t**3)
    )
    # The iteration runs 20 digits above the caller's precision, so the stopping
    # test 10^-(dps + 5) is reachable and the result is good to the caller's.
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    with mpmath.extradps(20):
        log_p = mpmath.log(p)
        for _ in range(100):
            cdf = _mp_norm_cdf(x)
            step = (mpmath.log(cdf) - log_p) * cdf / mpmath.npdf(x)
            x -= step
            if abs(step) <= eps * abs(x):
                return x
    raise ArithmeticError(f"norm_ppf({p}) did not converge")


NORM_PPF = Spec(
    name="norm_ppf",
    description="Standard normal quantile function Phi^-1(p) (mean 0, sd 1)",
    inputs=("p",),
    tolerance=1e-13,
    tolerance_reason=(
        "The input p is an exact double and the quantities a quantile method needs, p - 1/2 "
        "(exact for 1/4 <= p <= 1 by Sterbenz's lemma), 1 - p (exact for p >= 1/2) and "
        "log p, are exact or correctly rounded, so nothing beyond a few ulps of arithmetic is "
        "unavoidable (Wichura's AS 241 reaches about 1e-16). The quantile is "
        "ill-conditioned in p near p = 1/2 (|p / (x phi(x))| grows without bound) but that "
        "only matters for an inexact p. 1e-13 as for the other functions; SciPy's worst on "
        "this grid is 2.1e-16"
    ),
    grid_description=(
        "p = 10^-k for k = 1..20 and 30..300 step 10, the doubles nearest 1 - 10^-k for "
        "k = 1..16, i/64 for i = 0..64, 1/2 +- 10^-k for k = 5, 10, 15, and 5e-324, 1e-310, "
        "2^-1022; p = 0 and p = 1 are tagged pole (the quantile is -inf and +inf)"
    ),
    grid=_norm_ppf_grid,
    mpmath_fn=_mp_norm_ppf,
    scipy_name="scipy.stats.norm.ppf",
    scipy_fn=st.norm.ppf,
    input_tags=_norm_ppf_tags,
)


# --- Student's t distribution ------------------------------------------------

# Degrees of freedom: heavy tails (0.1 .. 2.5), small integers, and large values
# where t approaches the standard normal.
_T_DF = (0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 50.0, 100.0,
         1e3, 1e4, 1e5, 1e6, 1e8)  # fmt: skip
_T_X = (1e-300, 1e-10, 0.01, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,
        15.0, 20.0, 30.0, 40.0, 100.0, 1e3, 1e5, 1e10, 1e20, 1e50, 1e100, 1e300)  # fmt: skip


def _t_x_grid() -> list[tuple[float, float]]:
    xs = sorted({0.0, *_T_X, *(-x for x in _T_X)})
    return [(df, x) for df in _T_DF for x in xs]


def _mp_betainc_cf(a: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """I_x(a, b) for 0 < x < 1 by the continued fraction DLMF 8.17.22,

        I_x(a, b) = x^a (1-x)^b / (a B(a, b)) * 1 / (1 + d_1 / (1 + d_2 / (1 + ...))),

    d_(2m+1) = -(a+m)(a+b+m) x / ((a+2m)(a+2m+1)), d_(2m) = m(b-m) x / ((a+2m-1)(a+2m)),
    evaluated with the modified Lentz method for x <= (a+1)/(a+b+2), where it
    converges fast, and through I_x(a, b) = 1 - I_(1-x)(b, a) otherwise.

    The t datasets use it instead of _mp_betainc, whose series needs about
    1/(1-x) terms (1e8 for df = 1e8 at |t| = 1), and instead of mpmath.betainc,
    whose hypergeometric evaluation fails to converge for df >= 1e5. Tested
    against both and against closed forms in tests/test_functions.py.
    """
    y = mpmath.fsub(1, x, exact=True)
    if x > (a + 1) / (a + b + 2):
        return 1 - _mp_betainc_cf(b, a, y)
    eps = mpmath.mpf(2) ** (-mpmath.mp.prec - 10)
    tiny = mpmath.mpf(2) ** (-10 * mpmath.mp.prec)
    with mpmath.extraprec(30):
        f = c = mpmath.mpf(1)
        d = mpmath.mpf(0)
        for n in range(1, 10**7):
            m = n // 2
            if n % 2:
                term = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
            else:
                term = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
            d = 1 + term * d
            d = 1 / (d if d != 0 else tiny)
            c = 1 + term / c
            c = c if c != 0 else tiny
            delta = c * d
            f *= delta
            if abs(delta - 1) < eps:
                break
        else:
            raise ArithmeticError(f"betainc continued fraction ({a}, {b}, {x}) did not converge")
        # f converges to 1 + d_1 / (1 + d_2 / (1 + ...)) (Lentz from b_0 = 1).
        return x**a * y**b / (a * mpmath.beta(a, b)) / f


def _mp_t_lower_tail(df: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(-|x|) for the t distribution with df degrees of freedom (A&S 26.7.1,
    DLMF 8.17.2): F(-|x|) = I_h(df/2, 1/2) / 2 with h = df / (df + x^2).

    h is exact at the working precision for exact double inputs except when
    x^2 / df is below the precision (then F rounds to 1/2 as a double anyway);
    z = 1 - h = x^2 / (df + x^2) is passed separately, so the complement
    I_h(a, b) = 1 - I_z(b, a), taken where h is above the mean, never forms
    1 minus a number near 1.
    """
    if x == 0:
        return mpmath.mpf(1) / 2
    half = mpmath.mpf(1) / 2
    a = df / 2
    x2 = x * x
    h = df / (df + x2)
    z = x2 / (df + x2)
    if h > (a + 1) / (a + half + 2):
        return (1 - _mp_betainc_cf(half, a, z)) / 2
    return _mp_betainc_cf(a, half, h) / 2


def _mp_t_cdf(df: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    tail = _mp_t_lower_tail(df, x)
    return tail if x <= 0 else 1 - tail


def _mp_t_sf(df: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    return _mp_t_cdf(df, -x)


def _mp_t_pdf(df: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # A&S 26.7.1.
    return (1 + x * x / df) ** (-(df + 1) / 2) / (mpmath.sqrt(df) * mpmath.beta(df / 2, 0.5))


def t_cdf_condition(df: float, x: float) -> mpmath.mpf:
    """|x f(x) / F(x)|: the relative change of F per relative change of x.

    About x^2 where t is close to normal and tends to df in the power-law tail
    (F(x) is about C |x|^-df there); the survival function has the same
    condition number at -x.
    """
    df, x = mpmath.mpf(df), mpmath.mpf(x)
    return abs(x * _mp_t_pdf(df, x) / _mp_t_cdf(df, x))


def _t_ppf_grid() -> list[tuple[float, float]]:
    # Both tails (p = 10^-k and the doubles nearest 1 - 10^-k), the smallest
    # subnormal, the bulk, and points just off 1/2; p = 0 and 1 are poles.
    lower = [5e-324, 1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4]
    upper = [0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15]
    middle = [0.5 - 1e-10, 0.5, 0.5 + 1e-10]
    ps = sorted({0.0, 1.0, *lower, *middle, *upper})
    return [(df, p) for df in _T_DF for p in ps]


def _t_ppf_tags(df: float, p: float) -> list[str]:
    # The quantile is -inf at p = 0 and +inf at p = 1.
    return [POLE] if p in (0.0, 1.0) else []


def _mp_t_ppf(df: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with F(x) = p (0 < p < 1), by Newton's method on log F(-e^u).

    p is an exact double, so 1 - p is exact at any working precision above 53
    bits and the upper half uses the symmetry ppf(p) = -ppf(1 - p). In u =
    log|x| the function log F(-e^u) is smooth and decreasing, with derivative
    -|x| f(x) / F(x) (the condition number); bisection brackets the root to
    within 1e-3 and Newton's method then converges quadratically.
    """
    half = mpmath.mpf(1) / 2
    if p == half:
        return mpmath.mpf(0)
    if p > half:
        return -_mp_t_ppf(df, 1 - p)
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    with mpmath.extradps(20):
        log_p = mpmath.log(p)

        def g(u: mpmath.mpf) -> mpmath.mpf:
            return mpmath.log(_mp_t_lower_tail(df, mpmath.exp(u))) - log_p

        lo, hi = mpmath.mpf(-1), mpmath.mpf(1)
        while g(lo) <= 0:
            lo *= 2
        while g(hi) >= 0:
            hi *= 2
        while hi - lo > mpmath.mpf(1) / 1000:
            mid = (lo + hi) / 2
            if g(mid) > 0:
                lo = mid
            else:
                hi = mid
        u = (lo + hi) / 2
        for _ in range(100):
            x = mpmath.exp(u)
            slope = -x * _mp_t_pdf(df, x) / _mp_t_lower_tail(df, x)
            step = g(u) / slope
            u -= step
            if abs(step) <= eps:
                return -mpmath.exp(u)
    raise ArithmeticError(f"t_ppf({df}, {p}) did not converge")


_T_DF_DESCRIPTION = (
    "df in {0.1, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 20, 50, 100, 1e3, 1e4, 1e5, 1e6, 1e8}"
)

T_CDF = Spec(
    name="t_cdf",
    description="Student's t cumulative distribution function F(x; df) (location 0, scale 1)",
    inputs=("df", "x"),
    tolerance=2e-13,
    tolerance_reason=(
        "Rounding a quantity proportional to x or x^2 (such as x^2 / df) to a double changes F "
        "by the condition number |x f(x) / F(x)| times 2^-53. It is about x^2 where t is close "
        "to normal and tends to df in the power-law tail; its largest value at a normal result "
        "on this grid is 901 (df = 1e8, x = -30), so up to 1.0e-13 is unavoidable for a method "
        "that rounds such a quantity; 2e-13 allows about 2x that, as for norm_cdf. SciPy's worst "
        "error here is 1.03e-13 (df = 1e4, x = -30), except that it returns 0.5 at x = +-1e-10 "
        "for df = 1 (error 6.4e-11) and 0 at x = -1e300 for df <= 1. Results below 2.2e-308 "
        "are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_T_DF_DESCRIPTION} x x = 0 and +-{{1e-300, 1e-10, 0.01, 0.1, 0.25, 0.5, 1, 1.5, 2, "
        "3, 4, 5, 7.5, 10, 15, 20, 30, 40, 100, 1e3, 1e5, 1e10, 1e20, 1e50, 1e100, 1e300}; the "
        "lower tail decays like |x|^-df, so it is 0 as a double for large |x| unless df is small"
    ),
    grid=_t_x_grid,
    mpmath_fn=_mp_t_cdf,
    scipy_name="scipy.stats.t.cdf",
    scipy_fn=lambda df, x: st.t.cdf(x, df),
)

T_SF = Spec(
    name="t_sf",
    description="Student's t survival function 1 - F(x; df) = F(-x; df) (location 0, scale 1)",
    inputs=("df", "x"),
    tolerance=2e-13,
    tolerance_reason=(
        "1 - F(x) = F(-x) and the grid is symmetric, so the reasoning of t_cdf applies "
        "mirrored: the condition number |x f(x) / (1 - F(x))| reaches 901 at a normal result "
        "(df = 1e8, x = 30), so up to 1.0e-13 is unavoidable; 2e-13 allows about 2x that. SciPy's "
        "worst error here is 1.03e-13 (df = 1e4, x = 30), except that it returns 0.5 at "
        "x = +-1e-10 for df = 1 and 0 at x = 1e300 for df <= 1"
    ),
    grid_description=f"the t_cdf grid ({_T_DF_DESCRIPTION} x 53 values of x)",
    grid=_t_x_grid,
    mpmath_fn=_mp_t_sf,
    scipy_name="scipy.stats.t.sf",
    scipy_fn=lambda df, x: st.t.sf(x, df),
)

T_PPF = Spec(
    name="t_ppf",
    description="Student's t quantile function F^-1(p; df) (location 0, scale 1)",
    inputs=("df", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "p and df are exact doubles and 1 - p is exact for p >= 1/2, but in the tail the "
        "quantile is |x| = exp(L) with L = log|x| of the order of log(1/p) / df; rounding L "
        "(or 1/df, which multiplies it) to a double changes x by |L| times 2^-53. The largest "
        "|L| at a finite, normal quantile on this grid is 689.6 (df = 1, p = 1e-300, "
        "x = -3.2e299), so up to 7.7e-14 is unavoidable; 2e-13 allows 2.6x that and matches "
        "t_cdf. SciPy's worst error at the points it gets within tolerance is 1.3e-14; it fails "
        "36 points: up to 4.3e-7 off at p = 1/2 +- 1e-10 for 1 <= df <= 20 (0 for df = 4), "
        "5e-4 to 9e-4 off at p = 5e-324 for df >= 1000, +inf (wrong sign) at p = 5e-324 for "
        "2.5 <= df <= 100 (also at 1e-300 for df <= 10), and quantiles capped near -1e153 for "
        "df <= 1.5"
    ),
    grid_description=(
        f"{_T_DF_DESCRIPTION} x p in {{0, 5e-324, 1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, "
        "1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, "
        "0.99, 0.999, and the doubles nearest 1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1}; p = 0 and 1 "
        "are tagged pole (the quantile is -inf and +inf); quantiles beyond the largest double "
        "(df <= 1 with tiny p) are tagged overflow"
    ),
    grid=_t_ppf_grid,
    mpmath_fn=_mp_t_ppf,
    scipy_name="scipy.stats.t.ppf",
    scipy_fn=lambda df, p: st.t.ppf(p, df),
    input_tags=_t_ppf_tags,
)

# --- gamma distribution (shape a, rate 1) -------------------------------------

# Shapes: below 1 (density infinite at 0), small integers and half-integers
# (chi-squared with df = 2a), and large a where the distribution is close to
# normal with mean a and sd sqrt(a).
_GAMMA_A = (0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0,
            1e3, 1e4, 1e5, 1e6)  # fmt: skip
_GAMMA_X = (0.0, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0,
            5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 700.0, 750.0, 1e3, 1e4, 1e5,
            1e6, 1e7, 1e10, 1e300)  # fmt: skip
# Points a + k sqrt(a) (rounded to doubles) across the bulk and both tails.
_GAMMA_K = (-30, -10, -3, -1, 1, 3, 10, 30)


def _gamma_x_grid() -> list[tuple[float, float]]:
    points = set()
    for a in _GAMMA_A:
        points.update((a, x) for x in _GAMMA_X)
        points.update((a, x) for k in _GAMMA_K if (x := a + k * math.sqrt(a)) > 0)
    return sorted(points)


def _mp_gamma_prefix(a: mpmath.mpf, x: mpmath.mpf, shift: int) -> mpmath.mpf:
    # x^a e^-x / Gamma(a + shift), formed from its logarithm. The terms of the
    # logarithm reach 1.4e7 (a = 1e6) and cancel to about -k^2 / 2, so the
    # caller's 60 extra bits keep the result far beyond the working precision.
    return mpmath.exp(a * mpmath.log(x) - x - mpmath.loggamma(a + shift))


def _mp_gamma_lower_series(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """P(a, x) by the positive-term series DLMF 8.7.1 / 8.11.4,

        P(a, x) = x^a e^-x / Gamma(a + 1) * sum_k x^k / ((a + 1) ... (a + k)).

    Every term is positive, so there is no cancellation; the terms shrink once
    a + k > x. Used for x < a + 1, where about sqrt(a) * 15 terms suffice.
    """
    eps = mpmath.mpf(2) ** (-mpmath.mp.prec - 10)
    with mpmath.extraprec(60):
        term = total = mpmath.mpf(1)
        n = 0
        while True:
            n += 1
            term *= x / (a + n)
            total += term
            if a + n > x and term < eps * total:
                return _mp_gamma_prefix(a, x, 1) * total


def _mp_gamma_upper_cf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Q(a, x) by the continued fraction DLMF 8.9.2 (Legendre), in the even form

        Q(a, x) = x^a e^-x / Gamma(a) * 1 / (x + 1 - a - 1 (1 - a) / (x + 3 - a - 2 (2 - a) / ...)),

    evaluated with the modified Lentz method. Used for x >= a + 1, where it
    converges fast (about sqrt(a) * 10 steps near the mean).
    """
    eps = mpmath.mpf(2) ** (-mpmath.mp.prec - 10)
    tiny = mpmath.mpf(2) ** (-10 * mpmath.mp.prec)
    with mpmath.extraprec(60):
        b = x + 1 - a
        c = 1 / tiny
        d = 1 / b
        f = d
        for n in range(1, 10**7):
            an = -n * (n - a)
            b += 2
            d = an * d + b
            d = 1 / (d if d != 0 else tiny)
            c = b + an / c
            c = c if c != 0 else tiny
            delta = c * d
            f *= delta
            if abs(delta - 1) < eps:
                return _mp_gamma_prefix(a, x, 0) * f
    raise ArithmeticError(f"gamma continued fraction ({a}, {x}) did not converge")


def _mp_gamma_cdf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """P(a, x), the gamma distribution's cdf at x (shape a, rate 1).

    mpmath.gammainc fails to converge for a = 1e6 near the mean ("series
    converges too slowly"), so the dataset uses its own two routes, tested
    against mpmath.gammainc, finite sums and closed forms. The complements
    1 - Q and 1 - P are taken only where they are not small (x >= a + 1 and
    x < a + 1), so no digits are lost beyond the 60 extra bits.
    """
    if x == 0:
        return mpmath.mpf(0)
    if x < a + 1:
        return _mp_gamma_lower_series(a, x)
    return 1 - _mp_gamma_upper_cf(a, x)


def _mp_gamma_sf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Q(a, x) = 1 - P(a, x), the survival function (shape a, rate 1)."""
    if x == 0:
        return mpmath.mpf(1)
    if x < a + 1:
        return 1 - _mp_gamma_lower_series(a, x)
    return _mp_gamma_upper_cf(a, x)


def _mp_gamma_pdf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    return mpmath.exp((a - 1) * mpmath.log(x) - x - mpmath.loggamma(a))


def gamma_cdf_condition(a: float, x: float) -> mpmath.mpf:
    """|x f(x) / P(a, x)|: the relative change of P per relative change of x.

    About a in the lower tail (P is about x^a / Gamma(a + 1)) and about
    k sqrt(a) at x = a - k sqrt(a) for large a.
    """
    a, x = mpmath.mpf(a), mpmath.mpf(x)
    return abs(x * _mp_gamma_pdf(a, x) / _mp_gamma_cdf(a, x))


def gamma_sf_condition(a: float, x: float) -> mpmath.mpf:
    """|x f(x) / Q(a, x)|: about x in the upper tail (Q is about x^(a-1) e^-x / Gamma(a))."""
    a, x = mpmath.mpf(a), mpmath.mpf(x)
    return abs(x * _mp_gamma_pdf(a, x) / _mp_gamma_sf(a, x))


def _gamma_ppf_grid() -> list[tuple[float, float]]:
    # The t_ppf probabilities: both tails, the smallest subnormal, the bulk
    # and 1/2 +- 1e-10; p = 0 and 1 are poles.
    lower = [5e-324, 1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4]
    upper = [0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15]
    middle = [0.5 - 1e-10, 0.5, 0.5 + 1e-10]
    ps = sorted({0.0, 1.0, *lower, *middle, *upper})
    return [(a, p) for a in _GAMMA_A for p in ps]


def _gamma_ppf_tags(a: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is 0, the end of the support:
    # a finite value, so not a pole (a crate returning 0 there is right).
    return [POLE] if p == 1.0 else []


def _mp_gamma_ppf(a: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with P(a, x) = p (0 <= p < 1): _mp_gamma_quantile on the lower tail
    for p <= 1/2 and on the upper tail, 1 - p (exact for a double p at any
    working precision above 53 bits), above. ppf(0) = 0."""
    if p == 0:
        return mpmath.mpf(0)
    lower = p <= mpmath.mpf(1) / 2
    with mpmath.extradps(20):
        z = -mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * p)  # the normal quantile
    return _mp_gamma_quantile(a, p if lower else 1 - p, lower, z)


def _mp_gamma_quantile(a: mpmath.mpf, q: mpmath.mpf, lower: bool, z: mpmath.mpf) -> mpmath.mpf:
    """The x with P(a, x) = q (lower) or Q(a, x) = q (upper), 0 < q <= 1/2, by
    Newton's method on a logarithm.

    The root of log P(a, e^u) = log q, or of log Q(a, e^u) = log q. In u = log x
    both sides are smooth and monotone, with slopes x f / P and -x f / Q. The
    search starts from the Wilson-Hilferty approximation with z, the normal
    quantile of P(a, x) (or, where that is not positive, from the lower-tail
    power law P = x^a / Gamma(a + 1); where z is +inf because P(a, x) = 1 - q
    rounds to 1, from x = a - log q, about where Q = q), brackets the root,
    bisects to within 1e-3 and then converges quadratically.
    """
    sign = 1 if lower else -1
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    with mpmath.extradps(20):
        log_q = mpmath.log(q)

        def g(u: mpmath.mpf) -> mpmath.mpf:
            # Increasing in u for both branches after the sign.
            x = mpmath.exp(u)
            tail = _mp_gamma_cdf(a, x) if lower else _mp_gamma_sf(a, x)
            return sign * (mpmath.log(tail) - log_q)

        w = 1 - 1 / (9 * a) + z / (3 * mpmath.sqrt(a))
        if w == mpmath.inf:
            u0 = mpmath.log(a - log_q)
        elif w > 0:
            u0 = 3 * mpmath.log(w) + mpmath.log(a)
        else:
            u0 = (log_q + mpmath.loggamma(a + 1)) / a
        lo = hi = u0
        step = mpmath.mpf(1) / 1000
        while g(lo) > 0:
            lo -= step
            step *= 2
        step = mpmath.mpf(1) / 1000
        while g(hi) < 0:
            hi += step
            step *= 2
        while hi - lo > mpmath.mpf(1) / 1000:
            mid = (lo + hi) / 2
            if g(mid) > 0:
                hi = mid
            else:
                lo = mid
        u = (lo + hi) / 2
        for _ in range(100):
            x = mpmath.exp(u)
            tail = _mp_gamma_cdf(a, x) if lower else _mp_gamma_sf(a, x)
            step = g(u) / (x * _mp_gamma_pdf(a, x) / tail)
            u -= step
            if abs(step) <= eps:
                return mpmath.exp(u)
    raise ArithmeticError(f"gamma quantile ({a}, {q}, lower={lower}) did not converge")


_GAMMA_A_DESCRIPTION = "a in {0.01, 0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 20, 50, 100, 1e3, 1e4, 1e5, 1e6}"

GAMMA_CDF = Spec(
    name="gamma_cdf",
    description=(
        "Gamma distribution cumulative distribution function P(a, x) (shape a, rate 1); "
        "chi-squared with df = 2a at 2x"
    ),
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "With rate 1 this is P(a, x), the function of the gammainc dataset, so it gets "
        "gammainc's tolerance (accepted by the owner, to be revisited). The condition number "
        "|x f(x) / P| reaches 3.0e4 (a = 1e6, x = a - 30 sqrt(a)), so a method that rounds x / a "
        "loses up to 3.3e-12 there; one that forms x - a exactly (Sterbenz's lemma, "
        "a/2 <= x <= 2a) and expands in (x - a) / a, as Temme's uniform expansion does (DLMF "
        "8.12; Boost), loses only about k^2 2^-53 at x = a + k sqrt(a) (1e-13 at k = 30), so "
        "1e-12 demands such a method. SciPy's worst error here is 3.5e-13 (a = 1000, x = 500), "
        "except 2.1e-10 at a = 1e6, x = 990000 (P = 5.4e-24; the finite sum DLMF 8.4.10 "
        "confirms mpmath) and a subnormal result flushed to 0 at a = 50, x = 1e-5. Results "
        "below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_GAMMA_A_DESCRIPTION} x x in {{0, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, "
        "0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 700, 750, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10, "
        "1e300}, plus x = a + k sqrt(a) (rounded to a double) for k in {-30, -10, -3, -1, 1, 3, "
        "10, 30} where positive"
    ),
    grid=_gamma_x_grid,
    mpmath_fn=_mp_gamma_cdf,
    scipy_name="scipy.stats.gamma.cdf",
    scipy_fn=lambda a, x: st.gamma.cdf(x, a),
)

GAMMA_SF = Spec(
    name="gamma_sf",
    description=(
        "Gamma distribution survival function Q(a, x) = 1 - P(a, x) (shape a, rate 1); "
        "chi-squared with df = 2a at 2x"
    ),
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "Q = 1 - P, with gamma_cdf's tolerance for the same reasons: the condition number "
        "|x f(x) / Q| reaches 3.0e4 at a = 1e6, x = a + 30 sqrt(a), and is about x in the upper "
        "tail (up to about 745 at a normal result). SciPy's worst error here is 7.4e-13 "
        "(a = 1000, x = a + 30 sqrt(a)), except a subnormal result flushed to 0 at a = 5, "
        "x = 750"
    ),
    grid_description="the gamma_cdf grid",
    grid=_gamma_x_grid,
    mpmath_fn=_mp_gamma_sf,
    scipy_name="scipy.stats.gamma.sf",
    scipy_fn=lambda a, x: st.gamma.sf(x, a),
)

GAMMA_PPF = Spec(
    name="gamma_ppf",
    description=(
        "Gamma distribution quantile function, the x with P(a, x) = p (shape a, rate 1); "
        "half the chi-squared quantile with df = 2a"
    ),
    inputs=("a", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "p and a are exact doubles and 1 - p is exact for p >= 1/2, but in the lower tail the "
        "quantile is x = exp(L) with L = (log p + log Gamma(a + 1)) / a approximately; rounding L "
        "(or 1/a, which multiplies it) to a double changes x by |L| times 2^-53. The largest |L| "
        "at a nonzero, normal quantile on this grid is 691.3 (a = 0.01, p = 1e-3, x = 5.7e-301), "
        "so up to 7.7e-14 is unavoidable; 2e-13 allows 2.6x that and matches t_ppf. SciPy's "
        "worst error at the points it gets within tolerance is 1.1e-13 (a = 0.01, p = 0.01); it "
        "fails 6 points: 4.2e-12 to 4.7e-6 off at p = 5e-324 for a >= 100, and 4.7e-11 off at "
        "a = 1e6, p = 1e-10"
    ),
    grid_description=(
        f"{_GAMMA_A_DESCRIPTION} x p in {{0, 5e-324, 1e-300, 1e-100, 1e-50, 1e-20, 1e-10, "
        "1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, "
        "0.95, 0.99, 0.999, and the doubles nearest 1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1}; the "
        "quantile is 0 at p = 0 (the end of the support) and p = 1 is tagged pole (+inf); for "
        "small a the lower-tail quantiles are below the smallest double"
    ),
    grid=_gamma_ppf_grid,
    mpmath_fn=_mp_gamma_ppf,
    scipy_name="scipy.stats.gamma.ppf",
    scipy_fn=lambda a, p: st.gamma.ppf(p, a),
    input_tags=_gamma_ppf_tags,
)

# --- inverse gamma distribution (shape a, rate 1) -----------------------------

# X is inverse gamma with shape a and rate 1 exactly when T = 1 / X is gamma
# with shape a and rate 1, so F(x) = Q(a, 1/x), 1 - F(x) = P(a, 1/x) and the
# quantile is 1 / (the gamma quantile at 1 - p). Shapes up to 1e4: see the
# tolerance reason. a = 150 was added after an off-grid scan: from about
# a = 108 to 171, x^(-a-1) overflows a double at the bulk's x while Gamma(a)
# does not, so a density formed as x^(-a-1) e^(-1/x) / Gamma(a) is inf or NaN.
_INVGAMMA_A = tuple(sorted({a for a in _GAMMA_A if a <= 1e4} | {150.0}))


def _invgamma_x_grid() -> list[tuple[float, float]]:
    # The reciprocals (rounded to doubles) of the gamma_cdf grid's t, so that
    # 1/x visits the same places, plus 0 and two x where 1/x overflows a double.
    points = set()
    for a in _INVGAMMA_A:
        ts = [t for t in _GAMMA_X if t > 0]
        ts += [t for k in _GAMMA_K if (t := a + k * math.sqrt(a)) > 0]
        points.update((a, 1 / t) for t in ts)
        points.update((a, x) for x in (0.0, 5e-324, 1e-310))
    return sorted(points)


def _mp_reciprocal(x: mpmath.mpf) -> mpmath.mpf:
    # 1/x with 60 extra bits, so rounding it costs nothing at the working precision.
    with mpmath.extraprec(60):
        return 1 / x


def _mp_invgamma_cdf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = Q(a, 1/x) (shape a, rate 1), 0 for x <= 0."""
    return _mp_gamma_sf(a, _mp_reciprocal(x)) if x > 0 else mpmath.mpf(0)


def _mp_invgamma_sf(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = P(a, 1/x) (shape a, rate 1), 1 for x <= 0."""
    return _mp_gamma_cdf(a, _mp_reciprocal(x)) if x > 0 else mpmath.mpf(1)


def invgamma_cdf_condition(a: float, x: float) -> mpmath.mpf:
    """|x f(x) / F(x)|, which equals gamma_sf_condition at t = 1/x: F(x) = Q(a, 1/x)
    and d log F / d log x = -d log Q / d log t."""
    return gamma_sf_condition(a, _mp_reciprocal(mpmath.mpf(x)))


def invgamma_sf_condition(a: float, x: float) -> mpmath.mpf:
    """|x f(x) / (1 - F(x))| = gamma_cdf_condition at t = 1/x."""
    return gamma_cdf_condition(a, _mp_reciprocal(mpmath.mpf(x)))


def _invgamma_ppf_grid() -> list[tuple[float, float]]:
    ps = sorted({p for _, p in _gamma_ppf_grid()})
    return [(a, p) for a in _INVGAMMA_A for p in ps]


def _mp_invgamma_ppf(a: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with F(x) = Q(a, 1/x) = p (0 <= p < 1): x = 1/t, with t the gamma
    quantile on the upper tail Q(a, t) = p for p <= 1/2 and on the lower tail
    P(a, t) = 1 - p above (1 - p exact). ppf(0) = 0, the end of the support."""
    if p == 0:
        return mpmath.mpf(0)
    upper = p <= mpmath.mpf(1) / 2
    with mpmath.extradps(20):
        # The normal quantile of P(a, t) = 1 - p; +inf where 1 - 2p rounds to 1.
        z = mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * p)
    t = _mp_gamma_quantile(a, p if upper else 1 - p, not upper, z)
    return _mp_reciprocal(t)


_INVGAMMA_A_DESCRIPTION = "a in {0.01, 0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 20, 50, 100, 150, 1e3, 1e4}"
_INVGAMMA_TOLERANCE_REASON = (
    "F(x) = Q(a, 1/x) and 1 - F(x) = P(a, 1/x), the functions of gamma_sf and gamma_cdf, so "
    "they get gamma's tolerance 1e-12, and the condition number |x f(x) / F(x)| is gamma's at "
    "t = 1/x. Unlike gamma's x, t = 1/x is not an exact double: rounding it (as SciPy, Boost "
    "and statrs all do) costs the condition number times 2^-53. Shapes stop at 1e4, where the "
    "largest condition number at a normal result is 3.0e3 (t = a -+ 30 sqrt(a)) and rounding "
    "costs at most 3.4e-13; at gamma's a = 1e5 and 1e6 it would cost 1.1e-12 and 3.3e-12, "
    "beyond the tolerance, so a correct method that rounds 1/x would fail there. Results below "
    "2.2e-308 are judged relative to that number (conventions). SciPy's worst error here is "
    "7.5e-13 (cdf, a = 1e3, x = 1/(a + 30 sqrt(a))), except one subnormal result flushed to 0 "
    "in each function (cdf at a = 5, x = 1/750; sf at a = 50, x = 1/1e-5)"
)
_INVGAMMA_X_DESCRIPTION = (
    f"{_INVGAMMA_A_DESCRIPTION} x x = 1/t (rounded to a double) for t in {{1e-300, 1e-100, "
    "1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 700, 750, "
    "1e3, 1e4, 1e5, 1e6, 1e7, 1e10, 1e300} and t = a + k sqrt(a) (rounded to a double first) "
    "for k in {-30, -10, -3, -1, 1, 3, 10, 30} where positive, plus x = 0, 5e-324 and 1e-310 "
    "(where 1/x overflows a double)"
)

INVGAMMA_CDF = Spec(
    name="invgamma_cdf",
    description=(
        "Inverse gamma distribution cumulative distribution function Q(a, 1/x) (shape a, rate 1)"
    ),
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=_INVGAMMA_TOLERANCE_REASON,
    grid_description=_INVGAMMA_X_DESCRIPTION,
    grid=_invgamma_x_grid,
    mpmath_fn=_mp_invgamma_cdf,
    scipy_name="scipy.stats.invgamma.cdf",
    scipy_fn=lambda a, x: st.invgamma.cdf(x, a),
)

INVGAMMA_SF = Spec(
    name="invgamma_sf",
    description=(
        "Inverse gamma distribution survival function P(a, 1/x) = 1 - F(x) (shape a, rate 1)"
    ),
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=_INVGAMMA_TOLERANCE_REASON,
    grid_description=f"the invgamma_cdf grid: {_INVGAMMA_X_DESCRIPTION}",
    grid=_invgamma_x_grid,
    mpmath_fn=_mp_invgamma_sf,
    scipy_name="scipy.stats.invgamma.sf",
    scipy_fn=lambda a, x: st.invgamma.sf(x, a),
)

INVGAMMA_PPF = Spec(
    name="invgamma_ppf",
    description=(
        "Inverse gamma distribution quantile function, the x with Q(a, 1/x) = p "
        "(shape a, rate 1); the reciprocal of the gamma quantile at 1 - p"
    ),
    inputs=("a", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "As for gamma_ppf: the quantile is x = 1/t with t = exp(L), and rounding L (or 1/a, which "
        "multiplies it in the gamma lower tail) to a double changes x by |log x| times 2^-53. "
        "The largest |log x| at a finite, normal quantile on this grid is 691.3 (a = 0.01, "
        "p = 0.999, x = 1.8e300), so up to 7.7e-14 is unavoidable; 2e-13 allows 2.6x that "
        "and matches gamma_ppf. Quantiles above the largest double are tagged overflow. SciPy "
        "fails 11 of the 15 points with p = 5e-324 (5.5e-13 to 1.6e-4 off; it is right at "
        "a = 1, 1.5, 2 and 3) and no other point"
    ),
    grid_description=(
        f"{_INVGAMMA_A_DESCRIPTION} x the gamma_ppf probabilities p in {{0, 5e-324, 1e-300, "
        "1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 1/2 - 1e-10, 1/2, "
        "1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, and the doubles nearest 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15, 1}; the quantile is 0 at p = 0 (the end of the support) and "
        "p = 1 is tagged pole (+inf); for small a the upper quantiles are above the largest "
        "double"
    ),
    grid=_invgamma_ppf_grid,
    mpmath_fn=_mp_invgamma_ppf,
    scipy_name="scipy.stats.invgamma.ppf",
    scipy_fn=lambda a, p: st.invgamma.ppf(p, a),
    input_tags=_gamma_ppf_tags,
)

# --- chi distribution (k degrees of freedom) -----------------------------------

# X is chi with k degrees of freedom exactly when X^2 / 2 is gamma with shape
# k/2 and rate 1, so F(x) = P(k/2, x^2/2), 1 - F(x) = Q(k/2, x^2/2) and the
# quantile is sqrt(2 t) with t the gamma quantile. k is an integer (statrs
# takes a u64). 160 and 161 straddle statrs's switch from the direct density
# to exp(ln_pdf); k stops at 1e4 (shape 5e3): see the tolerance reason.
_CHI_K = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 160.0, 161.0, 300.0, 1e3, 1e4)
# x where x * x / 2 underflows to 0 (5e-324, 1e-170), to one of the five
# smallest subnormals (4e-162: 1e-323) or to a subnormal (1e-160), or where
# x * x overflows (1e155, 1e300).
_CHI_X_EXTRA = (0.0, 5e-324, 1e-170, 4e-162, 1e-160, 1e155, 1e300)


def _chi_x_grid() -> list[tuple[float, float]]:
    # x = sqrt(2 t) (rounded to a double) for the gamma_cdf grid's t at shape
    # k/2, so that x^2/2 visits the same places, plus the extra x.
    points = set()
    for k in _CHI_K:
        a = k / 2
        ts = [t for t in _GAMMA_X if t > 0]
        ts += [t for j in _GAMMA_K if (t := a + j * math.sqrt(a)) > 0]
        points.update((k, math.sqrt(2 * t)) for t in ts)
        points.update((k, x) for x in _CHI_X_EXTRA)
    return sorted(points)


def _mp_half_square(x: mpmath.mpf) -> mpmath.mpf:
    # x^2 / 2 with 60 extra bits, so rounding it costs nothing at the working precision.
    with mpmath.extraprec(60):
        return x * x / 2


def _mp_chi_cdf(k: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = P(k/2, x^2/2), 0 for x <= 0."""
    return _mp_gamma_cdf(k / 2, _mp_half_square(x)) if x > 0 else mpmath.mpf(0)


def _mp_chi_sf(k: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = Q(k/2, x^2/2), 1 for x <= 0."""
    return _mp_gamma_sf(k / 2, _mp_half_square(x)) if x > 0 else mpmath.mpf(1)


def chi_cdf_condition(k: float, x: float) -> mpmath.mpf:
    """|x f(x) / F(x)| = 2 gamma_cdf_condition(k/2, x^2/2): F(x) = P(k/2, t) with
    t = x^2/2, and d log t / d log x = 2."""
    return 2 * gamma_cdf_condition(k / 2, _mp_half_square(mpmath.mpf(x)))


def chi_sf_condition(k: float, x: float) -> mpmath.mpf:
    """|x f(x) / (1 - F(x))| = 2 gamma_sf_condition(k/2, x^2/2)."""
    return 2 * gamma_sf_condition(k / 2, _mp_half_square(mpmath.mpf(x)))


def _chi_ppf_grid() -> list[tuple[float, float]]:
    ps = sorted({p for _, p in _gamma_ppf_grid()})
    return [(k, p) for k in _CHI_K for p in ps]


def _mp_chi_ppf(k: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with P(k/2, x^2/2) = p (0 <= p < 1): sqrt(2 t), t the gamma quantile."""
    t = _mp_gamma_ppf(k / 2, p)
    with mpmath.extraprec(60):
        return mpmath.sqrt(2 * t)


_CHI_K_DESCRIPTION = "k in {1, 2, 3, 5, 10, 20, 50, 100, 160, 161, 300, 1e3, 1e4}"
_CHI_TOLERANCE_REASON = (
    "F(x) = P(k/2, x^2/2) and 1 - F(x) = Q(k/2, x^2/2), the functions of gamma_cdf and "
    "gamma_sf, so they get gamma's tolerance 1e-12. The condition number |x f(x) / F(x)| is "
    "twice gamma's at t = x^2/2, but t is not an exact double: rounding x * x (as SciPy and "
    "statrs do; halving is exact) costs gamma's condition number at t times 2^-53. Degrees of "
    "freedom stop at 1e4 (shape 5e3), where gamma's largest condition number at a normal result "
    "is 2.1e3 (t = a -+ 30 sqrt(a)) and rounding costs at most 2.4e-13; at k = 2e6 (gamma's "
    "a = 1e6) it would cost 3.3e-12, beyond the tolerance. Results below 2.2e-308 are judged "
    "relative to that number (conventions). SciPy's worst error at the other points is 3.0e-13 "
    "(sf, k = 1000, x = sqrt(2 (a + 3 sqrt(a)))); it is 1.3e-12 and 5.1e-12 off at k = 1e4 "
    "(cdf at x = 75.88, sf at x = 119.34, where the finite sum DLMF 8.4.10 confirms mpmath), "
    "flushes 1 cdf and 3 sf subnormal results to 0, and for k = 1 loses digits where x * x "
    "is subnormal (x = 4e-162 and 1e-160) and returns 0 where it underflows (x = 1e-170)"
)
_CHI_X_DESCRIPTION = (
    f"{_CHI_K_DESCRIPTION} x x = sqrt(2 t) (rounded to a double) for t in {{1e-300, 1e-100, "
    "1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 700, 750, "
    "1e3, 1e4, 1e5, 1e6, 1e7, 1e10, 1e300} and t = a + j sqrt(a) with a = k/2 (rounded to a "
    "double first) for j in {-30, -10, -3, -1, 1, 3, 10, 30} where positive, plus x = 0, "
    "5e-324, 1e-170, 4e-162 and 1e-160 (where x * x / 2 underflows to 0, 1e-323 or a "
    "subnormal), 1e155 and 1e300 (where x * x overflows)"
)

CHI_CDF = Spec(
    name="chi_cdf",
    description=(
        "Chi distribution cumulative distribution function P(k/2, x^2/2) (k degrees of freedom)"
    ),
    inputs=("k", "x"),
    tolerance=1e-12,
    tolerance_reason=_CHI_TOLERANCE_REASON,
    grid_description=_CHI_X_DESCRIPTION,
    grid=_chi_x_grid,
    mpmath_fn=_mp_chi_cdf,
    scipy_name="scipy.stats.chi.cdf",
    scipy_fn=lambda k, x: st.chi.cdf(x, k),
)

CHI_SF = Spec(
    name="chi_sf",
    description=(
        "Chi distribution survival function Q(k/2, x^2/2) = 1 - F(x) (k degrees of freedom)"
    ),
    inputs=("k", "x"),
    tolerance=1e-12,
    tolerance_reason=_CHI_TOLERANCE_REASON,
    grid_description=f"the chi_cdf grid: {_CHI_X_DESCRIPTION}",
    grid=_chi_x_grid,
    mpmath_fn=_mp_chi_sf,
    scipy_name="scipy.stats.chi.sf",
    scipy_fn=lambda k, x: st.chi.sf(x, k),
)

CHI_PPF = Spec(
    name="chi_ppf",
    description=(
        "Chi distribution quantile function, the x with P(k/2, x^2/2) = p (k degrees of "
        "freedom); sqrt(2 t) with t the gamma quantile"
    ),
    inputs=("k", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "As for gamma_ppf: the quantile is x = sqrt(2 t) with t = exp(L), and rounding L (or 2/k, "
        "which multiplies it in the lower tail) to a double changes t by |L| times 2^-53 and x by "
        "half that. The largest |L| / 2 = |log(x / sqrt(2))| at a nonzero, normal quantile on "
        "this grid is 690.9 (k = 1, p = 1e-300, x = 1.25e-300), so up to 7.7e-14 is "
        "unavoidable; 2e-13 allows 2.6x that and matches gamma_ppf. Subnormal quantiles are "
        "judged relative to 2.2e-308 (conventions). SciPy's worst error at the points it gets "
        "within tolerance is 1.8e-14 (k = 3, p = 1e-300); it fails 4: 0 instead of 1.25e-300 "
        "at k = 1, p = 1e-300, and 1.2e-9 to 1.6e-8 off at p = 5e-324 for k >= 300"
    ),
    grid_description=(
        f"{_CHI_K_DESCRIPTION} x the gamma_ppf probabilities p in {{0, 5e-324, 1e-300, 1e-100, "
        "1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 1/2 - 1e-10, 1/2, "
        "1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, and the doubles nearest 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15, 1}; the quantile is 0 at p = 0 (the end of the support) and "
        "p = 1 is tagged pole (+inf)"
    ),
    grid=_chi_ppf_grid,
    mpmath_fn=_mp_chi_ppf,
    scipy_name="scipy.stats.chi.ppf",
    scipy_fn=lambda k, p: st.chi.ppf(p, k),
    input_tags=_gamma_ppf_tags,
)

# --- exponential distribution ---------------------------------------------------


def _expon_x_grid() -> list[tuple[float]]:
    # Tiny x (where 1 - e^-x cancels), x = 0..40 step 0.25 (exact; the cdf
    # rounds to 1 from about 37.4), the edge where e^-x becomes subnormal
    # (708.40) and 0 (745.13), decades for huge x, and negative x.
    small = [float(f"1e-{k}") for k in range(1, 21)] + [5e-324, 1e-310]
    dense = [i / 4 for i in range(161)]
    edge = [50.0, 100.0, 200.0, 500.0, 700.0, 705.0, 708.0, 708.25, 708.5, 709.0, 709.75,
            710.0, 720.0, 730.0, 740.0, 744.0, 745.0, 745.125, 745.25, 746.0, 750.0, 1000.0]  # fmt: skip
    decades = _decades()
    negative = [-0.25, -1.0, *(-d for d in decades)]
    return [(x,) for x in sorted(set(small + dense + edge + decades + negative))]


def _mp_expon_cdf(x: mpmath.mpf) -> mpmath.mpf:
    # F(x) = 1 - e^-x = -expm1(-x) for x >= 0 (rate 1), 0 below the support.
    # mpmath.expm1 keeps full relative precision for tiny x.
    return -mpmath.expm1(-x) if x > 0 else mpmath.mpf(0)


def _mp_expon_sf(x: mpmath.mpf) -> mpmath.mpf:
    return mpmath.exp(-x) if x > 0 else mpmath.mpf(1)


def expon_cdf_condition(x: float) -> mpmath.mpf:
    """|x f(x) / F(x)| = x / (e^x - 1): at most 1 (its limit at x = 0)."""
    x = mpmath.mpf(x)
    return x / mpmath.expm1(x)


def expon_sf_condition(x: float) -> mpmath.mpf:
    """|x f(x) / Q(x)| = x."""
    return abs(mpmath.mpf(x))


def _expon_ppf_tags(p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is 0, the end of the support:
    # a finite value, so not a pole.
    return [POLE] if p == 1.0 else []


def _mp_expon_ppf(p: mpmath.mpf) -> mpmath.mpf:
    # x = -log(1 - p) = -log1p(-p); -p is exact, and mpmath.log1p keeps full
    # relative precision for tiny p.
    return -mpmath.log1p(-p)


_EXPON_X_DESCRIPTION = (
    "x = 10^-k for k = 1..20, 5e-324, 1e-310, x = 0..40 step 0.25, 50, 100, 200, 500, 700, "
    "705, 708, 708.25, 708.5, 709, 709.75, 710, 720, 730, 740, 744, 745, 745.125, 745.25, 746, "
    "750, 1000, +-10^k for k = -300..300 step 10, -0.25 and -1"
)

EXPON_CDF = Spec(
    name="expon_cdf",
    description=(
        "Exponential distribution cumulative distribution function 1 - e^-x (rate 1); "
        "the gamma distribution with shape 1"
    ),
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "x is an exact double and -expm1(-x) of an exact argument is within an ulp or so in "
        "any good libm, so nothing beyond a few ulps is unavoidable. The condition number "
        "|x f(x) / F(x)| = x / (e^x - 1) is at most 1, so even a method that rounds a "
        "quantity proportional to x (rate times x for another rate) loses at most 2^-53. "
        "1e-13 as for the other functions; SciPy's worst error here is 2.2e-16 (one subnormal step "
        "at x = 5e-324 and 1e-310, judged relative to 2.2e-308)"
    ),
    grid_description=(
        f"{_EXPON_X_DESCRIPTION}; F is 0 for x <= 0 and rounds to 1 above about 37.4"
    ),
    grid=_expon_x_grid,
    mpmath_fn=_mp_expon_cdf,
    scipy_name="scipy.stats.expon.cdf",
    scipy_fn=st.expon.cdf,
)

EXPON_SF = Spec(
    name="expon_sf",
    description="Exponential distribution survival function e^-x (rate 1)",
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "x is an exact double and exp(-x) of an exact argument is within an ulp or so in any "
        "good libm. The condition number |x f(x) / Q(x)| = x reaches 708.4 at the last normal "
        "result, so a method that rounds a quantity proportional to x (rate times x for "
        "another rate, not checked here) would lose up to 7.9e-14 there; with rate 1 nothing "
        "is rounded. 1e-13 as for the other functions; SciPy's worst error here is 2.1e-16. "
        "Results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"the expon_cdf grid: {_EXPON_X_DESCRIPTION}; Q is 1 for x <= 0, subnormal from "
        "about x = 708.40 and 0 above about 745.13"
    ),
    grid=_expon_x_grid,
    mpmath_fn=_mp_expon_sf,
    scipy_name="scipy.stats.expon.sf",
    scipy_fn=st.expon.sf,
)

EXPON_PPF = Spec(
    name="expon_ppf",
    description="Exponential distribution quantile function -log(1 - p) (rate 1)",
    inputs=("p",),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double, so -p is exact and -log1p(-p) is within an ulp or so in any "
        "good libm (log(1 - p) also, for p >= 1/2, where 1 - p is exact); nothing beyond a "
        "few ulps is unavoidable. 1e-13 as for the other functions; SciPy's worst error here "
        "is 1.6e-16"
    ),
    grid_description=(
        "the norm_ppf grid: p = 10^-k for k = 1..20 and 30..300 step 10, the doubles nearest "
        "1 - 10^-k for k = 1..16, i/64 for i = 0..64, 1/2 +- 10^-k for k = 5, 10, 15, and "
        "5e-324, 1e-310, 2^-1022; the quantile is 0 at p = 0 (the end of the support) and "
        "p = 1 is tagged pole (+inf)"
    ),
    grid=_norm_ppf_grid,
    mpmath_fn=_mp_expon_ppf,
    scipy_name="scipy.stats.expon.ppf",
    scipy_fn=st.expon.ppf,
    input_tags=_expon_ppf_tags,
)

# --- Weibull distribution (shape k, scale lambda) ------------------------------

# F(x) = 1 - exp(-(x/lambda)^k) for x >= 0. Scales are powers of two, so
# x/lambda is exact for every grid x (each is u * lambda with u and the
# product normal doubles): only t = (x/lambda)^k needs rounding. Scales
# 2^-500 and 2^500 make lambda^-k overflow or underflow (to 0 or a
# subnormal) for shapes from 2 up, which a method that computes x^k and
# lambda^-k apart cannot survive.
_WEIBULL_SHAPES = (0.1, 0.5, 1.0, 1.5, 2.0, 2.1, 3.5, 10.0, 100.0)
_WEIBULL_SCALE_EXPONENTS = (-500, -10, 0, 10, 500)
# Targets for t = (x/lambda)^k: both tails, the bulk, where the cdf rounds to
# 1 (37), where the sf becomes subnormal (708.4) and 0 (745.13), and beyond.
_WEIBULL_T = (1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0,
              37.0, 50.0, 100.0, 300.0, 700.0, 708.5, 740.0, 745.5, 1e3, 1e10, 1e300)  # fmt: skip


def _weibull_x_grid() -> list[tuple[float, float, float]]:
    # x = u lambda with u the double nearest t^(1/k); kept where u and x are
    # normal and finite (then x/lambda = u exactly). Plus x = 0 and x = -lambda.
    points = set()
    for k in _WEIBULL_SHAPES:
        for m in _WEIBULL_SCALE_EXPONENTS:
            scale = math.ldexp(1.0, m)
            points.update({(k, scale, 0.0), (k, scale, -scale)})
            for t in _WEIBULL_T:
                # float() rounds to nearest; only normal u are kept.
                with mpmath.workdps(40):
                    u = float(mpmath.mpf(t) ** (1 / mpmath.mpf(k)))
                x = u * scale
                if 2.0**-1022 <= u < math.inf and 2.0**-1022 <= x < math.inf:
                    points.add((k, scale, x))
    return sorted(points)


def _mp_weibull_t(k: mpmath.mpf, scale: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # (x/lambda)^k with 60 extra bits, so rounding it costs nothing at the working precision.
    with mpmath.extraprec(60):
        return (x / scale) ** k


def _mp_weibull_cdf(k: mpmath.mpf, scale: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = -expm1(-(x/lambda)^k), 0 for x <= 0; mpmath.expm1 keeps full
    relative precision for tiny t."""
    return -mpmath.expm1(-_mp_weibull_t(k, scale, x)) if x > 0 else mpmath.mpf(0)


def _mp_weibull_sf(k: mpmath.mpf, scale: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = exp(-(x/lambda)^k), 1 for x <= 0."""
    return mpmath.exp(-_mp_weibull_t(k, scale, x)) if x > 0 else mpmath.mpf(1)


def weibull_cdf_condition_in_t(k: float, scale: float, x: float) -> mpmath.mpf:
    """|t dF/dt / F| = t / (e^t - 1) at t = (x/lambda)^k: at most 1."""
    t = _mp_weibull_t(mpmath.mpf(k), mpmath.mpf(scale), mpmath.mpf(x))
    return t / mpmath.expm1(t)


def weibull_sf_condition_in_t(k: float, scale: float, x: float) -> mpmath.mpf:
    """|t dQ/dt / Q| = t at t = (x/lambda)^k."""
    return _mp_weibull_t(mpmath.mpf(k), mpmath.mpf(scale), mpmath.mpf(x))


def _weibull_ppf_grid() -> list[tuple[float, float, float]]:
    ps = sorted({p for _, p in _gamma_ppf_grid()})
    return [
        (k, math.ldexp(1.0, m), p)
        for k in _WEIBULL_SHAPES
        for m in _WEIBULL_SCALE_EXPONENTS
        for p in ps
    ]


def _weibull_ppf_tags(k: float, scale: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is 0, the end of the support:
    # a finite value, so not a pole.
    return [POLE] if p == 1.0 else []


def _mp_weibull_ppf(k: mpmath.mpf, scale: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = lambda (-log1p(-p))^(1/k); -p is exact and mpmath.log1p keeps full
    relative precision for tiny p."""
    with mpmath.extraprec(60):
        return scale * (-mpmath.log1p(-p)) ** (1 / k)


_WEIBULL_PARAMETERS_DESCRIPTION = (
    "shape k in {0.1, 0.5, 1, 1.5, 2, 2.1, 3.5, 10, 100} x scale lambda = 2^m for m in "
    "{-500, -10, 0, 10, 500}"
)
_WEIBULL_X_DESCRIPTION = (
    f"{_WEIBULL_PARAMETERS_DESCRIPTION} x x = u lambda with u the double nearest t^(1/k) for "
    "t in {1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.1, 0.5, 1, 2, 5, 10, 20, 37, 50, 100, "
    "300, 700, 708.5, 740, 745.5, 1e3, 1e10, 1e300}, kept where u and x are normal and finite "
    "(so x/lambda = u exactly), plus x = 0 and x = -lambda"
)

WEIBULL_CDF = Spec(
    name="weibull_cdf",
    description=(
        "Weibull distribution cumulative distribution function 1 - exp(-(x/lambda)^k) "
        "(shape k, scale lambda)"
    ),
    inputs=("k", "lambda", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, k and lambda are exact doubles and lambda is a power of two, so x/lambda is exact "
        "and only t = (x/lambda)^k must be rounded (pow is within an ulp or so in any good "
        "libm), then -expm1(-t) is within an ulp. The condition number in t, t / (e^t - 1), is "
        "at most 1, so a careful method loses a few ulps at most. 1e-13 as for the other "
        "functions; results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_WEIBULL_X_DESCRIPTION}; F is 0 for x <= 0 and rounds to 1 above about t = 37.4"
    ),
    grid=_weibull_x_grid,
    mpmath_fn=_mp_weibull_cdf,
    scipy_name="scipy.stats.weibull_min.cdf",
    scipy_fn=lambda k, scale, x: st.weibull_min.cdf(x, k, scale=scale),
)

WEIBULL_SF = Spec(
    name="weibull_sf",
    description="Weibull distribution survival function exp(-(x/lambda)^k) (shape k, scale lambda)",
    inputs=("k", "lambda", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, k and lambda are exact doubles and lambda is a power of two, so x/lambda is exact "
        "and only t = (x/lambda)^k must be rounded. A good libm's pow is within about half "
        "an ulp (glibc's bound is 0.52 ulp, 1.04 2^-53 relative). The condition number in t "
        "is t itself, 708.4 at the last normal result, so that rounding costs up to 8.2e-14 "
        "there (7.9e-14 if correctly rounded, as for expon_sf with another rate); a pow "
        "error of a full ulp would cost 1.6e-13. SciPy's worst error here is 5.4e-14. "
        "Results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"the weibull_cdf grid: {_WEIBULL_X_DESCRIPTION}; Q is 1 for x <= 0, subnormal from "
        "about t = 708.40 and 0 above about t = 745.13"
    ),
    grid=_weibull_x_grid,
    mpmath_fn=_mp_weibull_sf,
    scipy_name="scipy.stats.weibull_min.sf",
    scipy_fn=lambda k, scale, x: st.weibull_min.sf(x, k, scale=scale),
)

WEIBULL_PPF = Spec(
    name="weibull_ppf",
    description="Weibull distribution quantile function lambda (-log(1 - p))^(1/k)",
    inputs=("k", "lambda", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double, so L = -log1p(-p) is within an ulp or so; raising it to 1/k "
        "multiplies that error by 1/k (at most 10 on this grid) and adds an ulp, and "
        "multiplying by lambda (a power of two) is exact. 1/k itself is not a double for most "
        "k: rounding it (or log(L) / k) costs |log(x/lambda)| 2^-53, at most 744.4 2^-53 = "
        "8.3e-14 at a normal result on this grid (p = 5e-324). 1e-13 as for the other "
        "functions; results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_WEIBULL_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities: p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1; the quantile is 0 at p = 0 (the end of the support) and "
        "p = 1 is tagged pole (+inf)"
    ),
    grid=_weibull_ppf_grid,
    mpmath_fn=_mp_weibull_ppf,
    scipy_name="scipy.stats.weibull_min.ppf",
    scipy_fn=lambda k, scale, p: st.weibull_min.ppf(p, k, scale=scale),
    input_tags=_weibull_ppf_tags,
)

# --- log-normal distribution (location mu, scale sigma) ------------------------

# F(x) = Phi((ln x - mu) / sigma) for x > 0. The grid puts w = (ln x - mu) /
# sigma on the norm_cdf steps, so both tails reach complete underflow; mu and
# sigma pairs cover a narrow, a shifted and a wide distribution, keeping
# |ln x| below 160 (x stays a normal double).
_LOGNORM_PARAMETERS = ((0.0, 0.25), (-3.0, 0.5), (0.0, 1.0), (2.0, 4.0))


def _lognorm_x_grid() -> list[tuple[float, float, float]]:
    # x = the double nearest exp(mu + sigma w) for w = -40..40 step 0.5, plus
    # x = 0 and x = -1 (below the support).
    points = set()
    for mu, sigma in _LOGNORM_PARAMETERS:
        points.update({(mu, sigma, 0.0), (mu, sigma, -1.0)})
        for i in range(-80, 81):
            with mpmath.workdps(40):
                x = float(mpmath.exp(mu + sigma * mpmath.mpf(i) / 2))
            points.add((mu, sigma, x))
    return sorted(points)


def _lognorm_w(mu: mpmath.mpf, sigma: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # (ln x - mu) / sigma with 60 extra bits, so rounding it costs nothing at the
    # working precision.
    with mpmath.extraprec(60):
        return (mpmath.log(x) - mu) / sigma


def _mp_lognorm_cdf(mu: mpmath.mpf, sigma: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = Phi((ln x - mu) / sigma), 0 for x <= 0."""
    return _mp_norm_cdf(_lognorm_w(mu, sigma, x)) if x > 0 else mpmath.mpf(0)


def _mp_lognorm_sf(mu: mpmath.mpf, sigma: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = Phi(-(ln x - mu) / sigma), 1 for x <= 0."""
    return _mp_norm_sf(_lognorm_w(mu, sigma, x)) if x > 0 else mpmath.mpf(1)


def lognorm_rounding_cost(mu: float, sigma: float, x: float, upper: bool) -> mpmath.mpf:
    """Relative change of F (upper: of 1 - F) when ln x and w = (ln x - mu) / sigma
    are each rounded to a double: the condition number in w, phi(w) / Phi(+-w),
    times |ln x| / sigma + |w|, in units of 2^-53."""
    mu, sigma, x = mpmath.mpf(mu), mpmath.mpf(sigma), mpmath.mpf(x)
    w = _lognorm_w(mu, sigma, x)
    tail = _mp_norm_sf(w) if upper else _mp_norm_cdf(w)
    return mpmath.npdf(w) / tail * (abs(mpmath.log(x)) / sigma + abs(w))


def _lognorm_ppf_grid() -> list[tuple[float, float, float]]:
    ps = [p for (p,) in _norm_ppf_grid()]
    return [(mu, sigma, p) for mu, sigma in _LOGNORM_PARAMETERS for p in ps]


def _lognorm_ppf_tags(mu: float, sigma: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is 0, the end of the support:
    # a finite value, so not a pole.
    return [POLE] if p == 1.0 else []


def _mp_lognorm_ppf(mu: mpmath.mpf, sigma: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = exp(mu + sigma Phi^-1(p)), 0 at p = 0."""
    if p == 0:
        return mpmath.mpf(0)
    with mpmath.extraprec(60):
        return mpmath.exp(mu + sigma * _mp_norm_ppf(p))


_LOGNORM_PARAMETERS_DESCRIPTION = "(mu, sigma) in {(0, 0.25), (-3, 0.5), (0, 1), (2, 4)}"
_LOGNORM_X_DESCRIPTION = (
    f"{_LOGNORM_PARAMETERS_DESCRIPTION} x x = the double nearest exp(mu + sigma w) for "
    "w = -40..40 step 0.5, plus x = 0 and x = -1"
)

LOGNORM_CDF = Spec(
    name="lognorm_cdf",
    description=(
        "Log-normal distribution cumulative distribution function Phi((ln x - mu) / sigma) "
        "(location mu, scale sigma)"
    ),
    inputs=("mu", "sigma", "x"),
    tolerance=8e-13,
    tolerance_reason=(
        "x is an exact double but ln x is not: rounding ln x costs |ln x| 2^-53 in w = "
        "(ln x - mu) / sigma after dividing by sigma, and rounding w costs |w| 2^-53. Phi's "
        "condition number in w, phi(w) / Phi(w), is about |w| in the lower tail, so a method "
        "that rounds ln x and w loses (|ln x| / sigma + |w|) phi(w) / Phi(w) 2^-53, at most "
        "3.4e-13 at a normal result on this grid (mu = -3, sigma = 0.5, w = -37.5); 8e-13 "
        "is 2.35x that, the margin norm_cdf keeps (4e-13 over its 1.6e-13: 2.5x) for "
        "roundings of up to an ulp rather than half an ulp. SciPy flushes the subnormal "
        "result at w = -38 to 0. Results below 2.2e-308 are judged relative to that number "
        "(conventions)"
    ),
    grid_description=(
        f"{_LOGNORM_X_DESCRIPTION}; F is 0 for x <= 0, subnormal from about w = -37.52, 0 "
        "below about w = -38.49 and rounds to 1 above about w = 8.3"
    ),
    grid=_lognorm_x_grid,
    mpmath_fn=_mp_lognorm_cdf,
    scipy_name="scipy.stats.lognorm.cdf",
    scipy_fn=lambda mu, sigma, x: st.lognorm.cdf(x, sigma, scale=math.exp(mu)),
)

LOGNORM_SF = Spec(
    name="lognorm_sf",
    description=(
        "Log-normal distribution survival function Phi(-(ln x - mu) / sigma) "
        "(location mu, scale sigma)"
    ),
    inputs=("mu", "sigma", "x"),
    tolerance=8e-13,
    tolerance_reason=(
        "The reasoning of lognorm_cdf mirrored: rounding ln x and w = (ln x - mu) / sigma "
        "costs (|ln x| / sigma + |w|) phi(w) / Phi(-w) 2^-53, about |w| (|ln x| / sigma + "
        "|w|) 2^-53 in the upper tail, at most 3.2e-13 at a normal result on this grid "
        "(mu = 2, sigma = 4, w = 37.5); 8e-13 as for lognorm_cdf (2.5x). SciPy flushes the "
        "subnormal result at w = 38 to 0. Results below 2.2e-308 are judged relative to that "
        "number (conventions)"
    ),
    grid_description=(
        f"the lognorm_cdf grid: {_LOGNORM_X_DESCRIPTION}; 1 - F is 1 for x <= 0, subnormal "
        "from about w = 37.52, 0 above about w = 38.49 and rounds to 1 below about w = -8.3"
    ),
    grid=_lognorm_x_grid,
    mpmath_fn=_mp_lognorm_sf,
    scipy_name="scipy.stats.lognorm.sf",
    scipy_fn=lambda mu, sigma, x: st.lognorm.sf(x, sigma, scale=math.exp(mu)),
)

LOGNORM_PPF = Spec(
    name="lognorm_ppf",
    description="Log-normal distribution quantile function exp(mu + sigma Phi^-1(p))",
    inputs=("mu", "sigma", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "The quantile is exp(y) with y = mu + sigma Phi^-1(p), so an absolute error in y is "
        "the relative error of the result. Rounding y costs |y| 2^-53, and a normal quantile "
        "q good to a few ulps (norm_ppf: within 2^-53 x 2.4 for statrs and SciPy) costs "
        "sigma |q| times that; |y| and sigma |q| are at most 152 and 154 on this grid "
        "(mu = 2, sigma = 4, p = 5e-324), so a careful method stays within about 6e-14. "
        "1e-13 as for the other quantiles; results below 2.2e-308 are judged relative to "
        "that number (conventions)"
    ),
    grid_description=(
        f"{_LOGNORM_PARAMETERS_DESCRIPTION} x the norm_ppf probabilities: p = 10^-k for "
        "k = 1..20 and 30..300 step 10, the doubles nearest 1 - 10^-k for k = 1..16, i/64 for "
        "i = 0..64, 1/2 +- 10^-k for k = 5, 10, 15, and 5e-324, 1e-310, 2^-1022; the "
        "quantile is 0 at p = 0 (the end of the support) and p = 1 is tagged pole (+inf)"
    ),
    grid=_lognorm_ppf_grid,
    mpmath_fn=_mp_lognorm_ppf,
    scipy_name="scipy.stats.lognorm.ppf",
    scipy_fn=lambda mu, sigma, p: st.lognorm.ppf(p, sigma, scale=math.exp(mu)),
    input_tags=_lognorm_ppf_tags,
)

# --- Pareto distribution (shape alpha, scale x_m) -------------------------------

# F(x) = 1 - (x_m/x)^alpha for x >= x_m. Scales are powers of two and every
# grid x is u x_m with u a double or a power of two times one, so x/x_m is
# exact as a real number (it need not be a double: u = 2^1100 is not).
# Shape 2^-5 makes 1/alpha = 32 exact, so its quantile carries no rounding of
# the exponent; its quantile at x_m = 2^-500 is moderate although
# (1 - p)^(-1/alpha) overflows.
_PARETO_SHAPES = (2.0**-5, 0.1, 0.5, 1.0, 1.5, 2.0, 3.5, 10.0, 100.0)
_PARETO_SCALE_EXPONENTS = (-500, -10, 0, 10, 500)
# u = x/x_m as (mantissa, power of two): just above the support start, where
# 1 - (x_m/x)^alpha cancels; the bulk; and far tails, where x_m/x is subnormal
# (1.1 2^-1030, 1.1 2^-1060) or below the smallest double (2^-1100, 2^-1500)
# while (x_m/x)^alpha is not for small shapes.
_PARETO_U = tuple((u, 0) for u in (
    1.0, 1 + 2.0**-52, 1 + 1e-15, 1 + 1e-10, 1 + 1e-5, 1.001, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0,
    1e5, 1e10, 1e20, 1e50, 1e100, 1e200, 1e300,
)) + ((1.1, 1030), (1.1, 1060), (1.0, 1100), (1.0, 1500))  # fmt: skip


def _pareto_x_grid() -> list[tuple[float, float, float]]:
    # x = ldexp(mantissa, e) x_m, kept where x is normal and finite (then x/x_m
    # is exactly mantissa 2^e). Below the support: x = x_m (1 - 2^-53), x_m/2,
    # 0 and -x_m.
    points = set()
    for alpha in _PARETO_SHAPES:
        for m in _PARETO_SCALE_EXPONENTS:
            xm = math.ldexp(1.0, m)
            below = (xm * (1 - 2.0**-53), xm / 2, 0.0, -xm)
            points.update((alpha, xm, x) for x in below)
            for mantissa, e in _PARETO_U:
                # frexp's exponent: mantissa 2^(m + e) is below 2^top (ldexp raises past 2^1024).
                top = math.frexp(mantissa)[1] + m + e
                x = math.ldexp(mantissa, m + e) if top <= 1024 else math.inf
                if 2.0**-1022 <= x < math.inf:
                    points.add((alpha, xm, x))
    return sorted(points)


def _mp_pareto_t(alpha: mpmath.mpf, xm: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # alpha log(x/x_m) with 60 extra bits; mpmath.log keeps full relative
    # precision for x/x_m just above 1 (x/x_m - 1 is exact at this precision).
    with mpmath.extraprec(60):
        return alpha * mpmath.log(x / xm)


def _mp_pareto_cdf(alpha: mpmath.mpf, xm: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = -expm1(-alpha log(x/x_m)), 0 for x < x_m; mpmath.expm1 keeps full
    relative precision for x just above x_m."""
    return -mpmath.expm1(-_mp_pareto_t(alpha, xm, x)) if x >= xm else mpmath.mpf(0)


def _mp_pareto_sf(alpha: mpmath.mpf, xm: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = (x_m/x)^alpha = exp(-alpha log(x/x_m)), 1 for x < x_m."""
    return mpmath.exp(-_mp_pareto_t(alpha, xm, x)) if x >= xm else mpmath.mpf(1)


def pareto_cdf_condition_in_t(alpha: float, xm: float, x: float) -> mpmath.mpf:
    """|t dF/dt / F| = t / (e^t - 1) at t = alpha log(x/x_m) > 0: at most 1."""
    t = _mp_pareto_t(mpmath.mpf(alpha), mpmath.mpf(xm), mpmath.mpf(x))
    return t / mpmath.expm1(t)


def pareto_sf_condition_in_t(alpha: float, xm: float, x: float) -> mpmath.mpf:
    """|t dQ/dt / Q| = t at t = alpha log(x/x_m)."""
    return _mp_pareto_t(mpmath.mpf(alpha), mpmath.mpf(xm), mpmath.mpf(x))


def _pareto_ppf_grid() -> list[tuple[float, float, float]]:
    ps = sorted({p for _, p in _gamma_ppf_grid()})
    return [
        (alpha, math.ldexp(1.0, m), p)
        for alpha in _PARETO_SHAPES
        for m in _PARETO_SCALE_EXPONENTS
        for p in ps
    ]


def _pareto_ppf_tags(alpha: float, xm: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is x_m, the start of the
    # support: a finite value, so not a pole.
    return [POLE] if p == 1.0 else []


def _mp_pareto_ppf(alpha: mpmath.mpf, xm: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = x_m exp(-log1p(-p) / alpha); -p is exact and mpmath.log1p keeps full
    relative precision for tiny p."""
    with mpmath.extraprec(60):
        return xm * mpmath.exp(-mpmath.log1p(-p) / alpha)


_PARETO_PARAMETERS_DESCRIPTION = (
    "shape alpha in {2^-5, 0.1, 0.5, 1, 1.5, 2, 3.5, 10, 100} x scale x_m = 2^m for m in "
    "{-500, -10, 0, 10, 500}"
)
_PARETO_X_DESCRIPTION = (
    f"{_PARETO_PARAMETERS_DESCRIPTION} x x = u x_m for u in {{1, 1 + 2^-52, 1 + 1e-15, "
    "1 + 1e-10, 1 + 1e-5, 1.001, 1.1, 1.5, 2, 5, 10, 100, 1e5, 1e10, 1e20, 1e50, 1e100, "
    "1e200, 1e300, 1.1 2^1030, 1.1 2^1060, 2^1100, 2^1500}, kept where x is normal and "
    "finite (so x/x_m = u exactly as a real number), plus x = x_m (1 - 2^-53), x_m/2, 0 "
    "and -x_m below the support"
)

PARETO_CDF = Spec(
    name="pareto_cdf",
    description=(
        "Pareto distribution cumulative distribution function 1 - (x_m/x)^alpha "
        "(shape alpha, scale x_m)"
    ),
    inputs=("alpha", "xm", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, alpha and x_m are exact doubles and x_m is a power of two, so x/x_m is exact "
        "(where it is a double; otherwise log x - log x_m is within an ulp of log(x/x_m) "
        "relative, as no cancellation occurs there) and x/x_m - 1 is exact for x/x_m <= 2; "
        "t = alpha log1p(x/x_m - 1) is then within a few ulps, and -expm1(-t) adds an ulp. "
        "The condition number in t, t / (e^t - 1), is at most 1, so a careful method loses "
        "a few ulps at most. 1e-13 as for the other functions; results below 2.2e-308 are "
        "judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_PARETO_X_DESCRIPTION}; F is 0 for x < x_m and rounds to 1 above about t = 37.4 "
        "with t = alpha log(x/x_m)"
    ),
    grid=_pareto_x_grid,
    mpmath_fn=_mp_pareto_cdf,
    scipy_name="scipy.stats.pareto.cdf",
    scipy_fn=lambda alpha, xm, x: st.pareto.cdf(x, alpha, scale=xm),
)

PARETO_SF = Spec(
    name="pareto_sf",
    description="Pareto distribution survival function (x_m/x)^alpha (shape alpha, scale x_m)",
    inputs=("alpha", "xm", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, alpha and x_m are exact doubles and x_m is a power of two. The condition number "
        "in t = alpha log(x/x_m) is t itself, 708.4 at the last normal result, so rounding "
        "t once costs up to 7.9e-14 there (as for weibull_sf); where x_m/x is a normal "
        "double, pow(x_m/x, alpha) costs only alpha ulps (at most 100, 1.1e-14) plus one. "
        "Results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"the pareto_cdf grid: {_PARETO_X_DESCRIPTION}; Q is 1 for x < x_m, subnormal from "
        "about t = 708.40 and 0 above about t = 745.13 with t = alpha log(x/x_m)"
    ),
    grid=_pareto_x_grid,
    mpmath_fn=_mp_pareto_sf,
    scipy_name="scipy.stats.pareto.sf",
    scipy_fn=lambda alpha, xm, x: st.pareto.sf(x, alpha, scale=xm),
)

PARETO_PPF = Spec(
    name="pareto_ppf",
    description="Pareto distribution quantile function x_m (1 - p)^(-1/alpha)",
    inputs=("alpha", "xm", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double; 1 - p is exact for p >= 1/2, and below that rounding it costs "
        "at most 2^-53 / alpha (3.6e-15 for alpha = 2^-5) in the result. x_m (1 - p)^(-1/alpha) "
        "is then one pow, within an ulp or so, and multiplying by x_m (a power of two) is "
        "exact. 1/alpha is exact for alpha in {2^-5, 0.5, 1, 2}; for the other shapes "
        "rounding it costs |log(x/x_m)| 2^-53, at most 345.4 2^-53 = 3.8e-14 at a normal "
        "result on this grid (alpha = 0.1, p = 1 - 1e-15). 1e-13 as for the other functions; "
        "results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_PARETO_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities: p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1; the quantile is x_m at p = 0 (the start of the "
        "support) and p = 1 is tagged pole (+inf)"
    ),
    grid=_pareto_ppf_grid,
    mpmath_fn=_mp_pareto_ppf,
    scipy_name="scipy.stats.pareto.ppf",
    scipy_fn=lambda alpha, xm, p: st.pareto.ppf(p, alpha, scale=xm),
    input_tags=_pareto_ppf_tags,
)

# --- Gumbel distribution (location mu, scale beta) ------------------------------

# F(x) = exp(-exp(-z)) with z = (x - mu) / beta. Scales are powers of two and
# every grid x is mu + beta z exactly with z a double, so x - mu = beta z is an
# exact double and z carries no rounding. The pairs cover the standard
# distribution, tiny and huge scales, two shifted ones and a location and scale
# near the largest double, where x - mu overflows although z is moderate.
_GUMBEL_PARAMETERS = (
    (0.0, 1.0), (0.0, 2.0**-500), (0.0, 2.0**500), (-10.0, 0.25), (100.0, 8.0),
    (-(2.0**1023), 2.0**1023),
)  # fmt: skip
# z: the bulk in steps of 1/4 (the lower tail, where F underflows below about
# z = -6.6, is steep); far tails (Q is subnormal from about z = 708.4 and 0
# above about z = 745.1); and points next to 0.
_GUMBEL_Z = tuple(
    sorted(
        {i / 4 for i in range(-28, 161)}
        | {-1e3, -100.0, -20.0, -10.0, -8.0, -7.5}
        | {50.0, 100.0, 200.0, 500.0, 700.0, 708.0, 708.5, 709.0, 710.0, 720.0, 740.0}
        | {745.0, 746.0, 800.0, 1e3, 1e5, 1e10, 1e100, 1e300}
        | {s * 2.0**-30 for s in (-1, 1)}
    )
)


def _gumbel_x_grid() -> list[tuple[float, float, float]]:
    # x = mu + beta z, kept where it is an exact, finite double that is 0 or
    # normal.
    points = []
    for mu, beta in _GUMBEL_PARAMETERS:
        for z in _GUMBEL_Z:
            # 2200 bits hold mu + beta z exactly for every double on this grid.
            with mpmath.workprec(2200):
                exact = mpmath.mpf(mu) + mpmath.mpf(beta) * mpmath.mpf(z)
                if abs(exact) >= mpmath.mpf(2) ** 1024:
                    continue
                x = float(exact)
                if mpmath.mpf(x) == exact and (x == 0 or abs(x) >= 2.0**-1022):
                    points.append((mu, beta, x))
    return sorted(points)


def _gumbel_z(mu: mpmath.mpf, beta: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # (x - mu) / beta with 60 extra bits (exact on the grid).
    with mpmath.extraprec(60):
        return (x - mu) / beta


def _mp_gumbel_cdf(mu: mpmath.mpf, beta: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = exp(-exp(-z)) with z = (x - mu) / beta."""
    with mpmath.extraprec(60):
        return mpmath.exp(-mpmath.exp(-_gumbel_z(mu, beta, x)))


def _mp_gumbel_sf(mu: mpmath.mpf, beta: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x) = -expm1(-exp(-z)); mpmath.expm1 keeps full relative precision in
    the upper tail, where exp(-z) is tiny."""
    with mpmath.extraprec(60):
        return -mpmath.expm1(-mpmath.exp(-_gumbel_z(mu, beta, x)))


def gumbel_cdf_condition_in_w(mu: float, beta: float, x: float) -> mpmath.mpf:
    """|w dF/dw / F| = w at w = exp(-z): F = exp(-w)."""
    return mpmath.exp(-_gumbel_z(mpmath.mpf(mu), mpmath.mpf(beta), mpmath.mpf(x)))


def gumbel_sf_condition_in_w(mu: float, beta: float, x: float) -> mpmath.mpf:
    """|w dQ/dw / Q| = w / (e^w - 1) at w = exp(-z): at most 1."""
    w = gumbel_cdf_condition_in_w(mu, beta, x)
    return w / mpmath.expm1(w)


# The quantile mu - beta log(-log p) is 0 at p = 1/e for mu = 0, where
# log(-log p) cancels: the doubles nearest 1/e and 1/e +- 10^-k probe it.
_GUMBEL_NEAR_INV_E = tuple(
    sorted(
        [float(mpmath.exp(-1))]
        + [float(mpmath.exp(-1) + s * mpmath.mpf(10) ** -k) for k in (5, 10, 15) for s in (-1, 1)]
    )
)


def _gumbel_ppf_grid() -> list[tuple[float, float, float]]:
    ps = sorted({p for _, p in _gamma_ppf_grid()} | set(_GUMBEL_NEAR_INV_E))
    return [(mu, beta, p) for mu, beta in _GUMBEL_PARAMETERS for p in ps]


def _gumbel_ppf_tags(mu: float, beta: float, p: float) -> list[str]:
    # The quantile is -inf at p = 0 and +inf at p = 1.
    return [POLE] if p in (0.0, 1.0) else []


def _mp_gumbel_ppf(mu: mpmath.mpf, beta: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = mu - beta log(-log p). Near p = 1/e, -log p is close to 1 and its log
    cancels: |1 + log p| is at least 1e-17 on this grid, so the cancellation
    costs about 57 of the 60 extra bits (tested)."""
    with mpmath.extraprec(60):
        return mu - beta * mpmath.log(-mpmath.log(p))


_GUMBEL_PARAMETERS_DESCRIPTION = (
    "(mu, beta) in {(0, 1), (0, 2^-500), (0, 2^500), (-10, 1/4), (100, 8), (-2^1023, 2^1023)}"
)
_GUMBEL_X_DESCRIPTION = (
    f"{_GUMBEL_PARAMETERS_DESCRIPTION} x x = mu + beta z for z = -7..40 step 1/4, z in "
    "{-1000, -100, -20, -10, -8, -7.5, 50, 100, 200, 500, 700, 708, 708.5, 709, 710, 720, "
    "740, 745, 746, 800, 1000, 1e5, 1e10, 1e100, 1e300} and z = +-2^-30, kept where x is an "
    "exact, finite double that is 0 or normal (so z = (x - mu) / beta is exact)"
)

GUMBEL_CDF = Spec(
    name="gumbel_cdf",
    description=(
        "Gumbel distribution cumulative distribution function exp(-exp(-(x - mu) / beta)) "
        "(location mu, scale beta)"
    ),
    inputs=("mu", "beta", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, mu and beta are exact doubles, beta is a power of two and x - mu = beta z with z "
        "a double, so z = (x - mu) / beta is exact. F = exp(-w) with w = exp(-z) has "
        "condition number w in w, at most 708.4 at a normal result, so rounding w once costs "
        "up to 7.9e-14 there (as for weibull_sf), and the outer exp adds an ulp. This assumes "
        "an exp within about half an ulp (as glibc's is); one good only to 1 ulp could lose "
        "about 1.5e-13 at the grid's worst point (w = 665), which would not by itself be a "
        "bug. 1e-13 as for the other functions; results below 2.2e-308 are judged relative "
        "to that number (conventions)"
    ),
    grid_description=(
        f"{_GUMBEL_X_DESCRIPTION}; F is subnormal from about z = -6.56, 0 below about "
        "z = -6.61 and rounds to 1 above about z = 36.7"
    ),
    grid=_gumbel_x_grid,
    mpmath_fn=_mp_gumbel_cdf,
    scipy_name="scipy.stats.gumbel_r.cdf",
    scipy_fn=lambda mu, beta, x: st.gumbel_r.cdf(x, loc=mu, scale=beta),
)

GUMBEL_SF = Spec(
    name="gumbel_sf",
    description=(
        "Gumbel distribution survival function -expm1(-exp(-(x - mu) / beta)) "
        "(location mu, scale beta)"
    ),
    inputs=("mu", "beta", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "z = (x - mu) / beta is exact on this grid (see gumbel_cdf). Q = -expm1(-w) with "
        "w = exp(-z) has condition number w / (e^w - 1), at most 1, in w, so rounding w and "
        "the expm1 cost a few ulps. 1e-13 as for the other functions; results below "
        "2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"the gumbel_cdf grid: {_GUMBEL_X_DESCRIPTION}; Q rounds to 1 below about z = -3.6, "
        "is subnormal from about z = 708.4 and 0 above about z = 745.1"
    ),
    grid=_gumbel_x_grid,
    mpmath_fn=_mp_gumbel_sf,
    scipy_name="scipy.stats.gumbel_r.sf",
    scipy_fn=lambda mu, beta, x: st.gumbel_r.sf(x, loc=mu, scale=beta),
)

GUMBEL_PPF = Spec(
    name="gumbel_ppf",
    description="Gumbel distribution quantile function mu - beta log(-log p)",
    inputs=("mu", "beta", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double. A careful method forms L = log(-log p) to a few ulps "
        "relative, also near p = 1/e where -log p is close to 1 (for example L = "
        "log1p(-log1p(d)) with d = e p - 1 formed in double-double); then mu - beta L costs "
        "(|mu| + |beta L|) / |mu - beta L| ulps, at most about 14 on this grid (mu = -10, "
        "beta = 1/4, p = 1 - 1e-15), 1.6e-15. 1e-13 as for the other quantiles"
    ),
    grid_description=(
        f"{_GUMBEL_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities (p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1) and the doubles nearest 1/e and 1/e +- 10^-k for "
        "k = 5, 10, 15, where the quantile is near 0 for mu = 0; p = 0 and 1 are tagged pole "
        "(-inf and +inf), quantiles beyond the largest double are tagged overflow"
    ),
    grid=_gumbel_ppf_grid,
    mpmath_fn=_mp_gumbel_ppf,
    scipy_name="scipy.stats.gumbel_r.ppf",
    scipy_fn=lambda mu, beta, p: st.gumbel_r.ppf(p, loc=mu, scale=beta),
    input_tags=_gumbel_ppf_tags,
)

# --- Laplace distribution (location mu, scale b) --------------------------------

# F(x) = e^z / 2 for z < 0 and 1 - e^-z / 2 for z >= 0, with z = (x - mu) / b.
# The Gumbel pairs, for the same reasons: scales are powers of two and every
# grid x is mu + b z exactly, so z is exact; the last pair puts x - mu above
# the largest double although z is moderate.
_LAPLACE_PARAMETERS = _GUMBEL_PARAMETERS
# z: the bulk in steps of 1/4 on both sides (F rounds to 1 above about
# z = 36.7, Q below about z = -36.7); both far tails (the smaller of F and Q is
# subnormal from about |z| = 707.7 and 0 from about |z| = 744.4); and points
# next to 0.
_LAPLACE_Z = tuple(
    sorted(
        {i / 4 for i in range(-160, 161)}
        | {
            s * t
            for s in (-1, 1)
            for t in (50.0, 100.0, 200.0, 500.0, 700.0, 707.0, 708.0, 709.0, 710.0, 720.0, 740.0,
                      744.0, 745.0, 746.0, 800.0, 1e3, 1e5, 1e10, 1e100, 1e300, 2.0**-30, 2.0**-60)
        }
    )
)  # fmt: skip


def _laplace_x_grid() -> list[tuple[float, float, float]]:
    # x = mu + b z, kept where it is an exact, finite double that is 0 or normal.
    points = []
    for mu, b in _LAPLACE_PARAMETERS:
        for z in _LAPLACE_Z:
            # 2200 bits hold mu + b z exactly for every double on this grid.
            with mpmath.workprec(2200):
                exact = mpmath.mpf(mu) + mpmath.mpf(b) * mpmath.mpf(z)
                if abs(exact) >= mpmath.mpf(2) ** 1024:
                    continue
                x = float(exact)
                if mpmath.mpf(x) == exact and (x == 0 or abs(x) >= 2.0**-1022):
                    points.append((mu, b, x))
    return sorted(points)


def _mp_laplace_cdf(mu: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = e^z / 2 for z < 0 and 1 - e^-z / 2 for z >= 0 (at least 1/2, so the
    subtraction does not cancel)."""
    with mpmath.extraprec(60):
        z = (x - mu) / b
        return mpmath.exp(z) / 2 if z < 0 else 1 - mpmath.exp(-z) / 2


def _mp_laplace_sf(mu: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Q(x) = F(mu - (x - mu)): the distribution is symmetric about mu."""
    with mpmath.extraprec(60):
        return _mp_laplace_cdf(mu, b, 2 * mu - x)


def _laplace_ppf_grid() -> list[tuple[float, float, float]]:
    # The gamma_ppf probabilities and the doubles next to 1/2, where the
    # quantile is mu + b log(2p) with 2p next to 1.
    ps = sorted({p for _, p in _gamma_ppf_grid()} | {0.5 - 2.0**-54, 0.5 + 2.0**-53})
    return [(mu, b, p) for mu, b in _LAPLACE_PARAMETERS for p in ps]


def _laplace_ppf_tags(mu: float, b: float, p: float) -> list[str]:
    # The quantile is -inf at p = 0 and +inf at p = 1.
    return [POLE] if p in (0.0, 1.0) else []


def _mp_laplace_ppf(mu: mpmath.mpf, b: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = mu + b log(2p) for p <= 1/2 and mu - b log(2 - 2p) for p > 1/2."""
    with mpmath.extraprec(60):
        return mu + b * mpmath.log(2 * p) if p <= 0.5 else mu - b * mpmath.log(2 - 2 * p)


def laplace_ppf_cancellation(mu: float, b: float, p: float) -> mpmath.mpf:
    """(|mu| + |b L|) / |mu + b L| with b L the quantile's log part: the ulps that
    adding mu to a rounded b L can cost."""
    m = mpmath.mpf(mu)
    x = _mp_laplace_ppf(m, mpmath.mpf(b), mpmath.mpf(p))
    return (abs(m) + abs(x - m)) / abs(x)


_LAPLACE_PARAMETERS_DESCRIPTION = (
    "(mu, b) in {(0, 1), (0, 2^-500), (0, 2^500), (-10, 1/4), (100, 8), (-2^1023, 2^1023)}"
)
_LAPLACE_X_DESCRIPTION = (
    f"{_LAPLACE_PARAMETERS_DESCRIPTION} x x = mu + b z for z = -40..40 step 1/4, "
    "z = +-{50, 100, 200, 500, 700, 707, 708, 709, 710, 720, 740, 744, 745, 746, 800, 1000, "
    "1e5, 1e10, 1e100, 1e300} and z = +-2^-30, +-2^-60, kept where x is an exact, finite double "
    "that is 0 or normal (so z = (x - mu) / b is exact)"
)

LAPLACE_CDF = Spec(
    name="laplace_cdf",
    description=(
        "Laplace distribution cumulative distribution function e^z / 2 (z < 0), "
        "1 - e^-z / 2 (z >= 0), z = (x - mu) / b (location mu, scale b)"
    ),
    inputs=("mu", "b", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, mu and b are exact doubles, b is a power of two and x - mu = b z with z a double, "
        "so z is exact; exp of an exact argument is within an ulp, and for z >= 0 the result "
        "is at least 1/2, so 1 - e^-z / 2 adds at most an ulp. 1e-13 as for the other "
        "functions; results below 2.2e-308 are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_LAPLACE_X_DESCRIPTION}; F is subnormal from about z = -707.7, 0 below about "
        "z = -744.4 and rounds to 1 above about z = 36.7"
    ),
    grid=_laplace_x_grid,
    mpmath_fn=_mp_laplace_cdf,
    scipy_name="scipy.stats.laplace.cdf",
    scipy_fn=lambda mu, b, x: st.laplace.cdf(x, loc=mu, scale=b),
)

LAPLACE_SF = Spec(
    name="laplace_sf",
    description=(
        "Laplace distribution survival function 1 - e^z / 2 (z < 0), e^-z / 2 (z >= 0), "
        "z = (x - mu) / b (location mu, scale b)"
    ),
    inputs=("mu", "b", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "the mirror image of laplace_cdf (Q(x) = F(2 mu - x)), with the same exact z and the "
        "same rounding. 1e-13 as for the other functions; results below 2.2e-308 are judged "
        "relative to that number (conventions)"
    ),
    grid_description=(
        f"the laplace_cdf grid: {_LAPLACE_X_DESCRIPTION}; Q rounds to 1 below about "
        "z = -36.7, is subnormal from about z = 707.7 and 0 above about z = 744.4"
    ),
    grid=_laplace_x_grid,
    mpmath_fn=_mp_laplace_sf,
    scipy_name="scipy.stats.laplace.sf",
    scipy_fn=lambda mu, b, x: st.laplace.sf(x, loc=mu, scale=b),
)

LAPLACE_PPF = Spec(
    name="laplace_ppf",
    description=(
        "Laplace distribution quantile function mu + b log(2p) (p <= 1/2), "
        "mu - b log(2 - 2p) (p > 1/2)"
    ),
    inputs=("mu", "b", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double, and so are 2p and 2 - 2p (Sterbenz), so a careful method forms "
        "L = log(2p) or log(2 - 2p) within an ulp; then mu + b L costs (|mu| + |b L|) / "
        "|mu + b L| ulps, at most about 14 on this grid (mu = 100, b = 8, p = 1e-5; tested), "
        "1.6e-15. 1e-13 as for the other quantiles"
    ),
    grid_description=(
        f"{_LAPLACE_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities (p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1) and the doubles next to 1/2 (1/2 - 2^-54, 1/2 + 2^-53); "
        "p = 0 and 1 are tagged pole (-inf and +inf), quantiles beyond the largest double are "
        "tagged overflow"
    ),
    grid=_laplace_ppf_grid,
    mpmath_fn=_mp_laplace_ppf,
    scipy_name="scipy.stats.laplace.ppf",
    scipy_fn=lambda mu, b, p: st.laplace.ppf(p, loc=mu, scale=b),
    input_tags=_laplace_ppf_tags,
)

# --- Cauchy distribution (location mu, scale b) ---------------------------------

# F(x) = 1/2 + atan(z) / pi with z = (x - mu) / b. The Gumbel pairs, for the
# same reasons: scales are powers of two and every grid x is mu + b z exactly,
# so z is exact; the last pair puts x - mu above the largest double although z
# is moderate.
_CAUCHY_PARAMETERS = _GUMBEL_PARAMETERS
# z: the bulk in steps of 1/4 on both sides; both heavy tails, where the smaller
# of F and Q is about 1 / (pi |z|) (so 1/2 + atan(z) / pi loses about log10 |z|
# digits, and all of them from about |z| = 1e16), down to subnormal from about
# |z| = 1.4e307; and points next to 0.
_CAUCHY_Z = tuple(
    sorted(
        {i / 4 for i in range(-160, 161)}
        | {
            s * t
            for s in (-1, 1)
            for t in (50.0, 100.0, 1e3, 1e4, 1e5, 1e8, 1e10, 1e15, 1e16, 1e17, 1e20, 1e50, 1e100,
                      1e300, 1e305, 1e308, 2.0**-30, 2.0**-60)
        }
    )
)  # fmt: skip


def _cauchy_x_grid() -> list[tuple[float, float, float]]:
    # x = mu + b z, kept where it is an exact, finite double that is 0 or normal.
    points = []
    for mu, b in _CAUCHY_PARAMETERS:
        for z in _CAUCHY_Z:
            # 2200 bits hold mu + b z exactly for every double on this grid.
            with mpmath.workprec(2200):
                exact = mpmath.mpf(mu) + mpmath.mpf(b) * mpmath.mpf(z)
                if abs(exact) >= mpmath.mpf(2) ** 1024:
                    continue
                x = float(exact)
                if mpmath.mpf(x) == exact and (x == 0 or abs(x) >= 2.0**-1022):
                    points.append((mu, b, x))
    return sorted(points)


def _mp_cauchy_cdf(mu: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = acot(-z) / pi for z < 0 and 1 - acot(z) / pi for z >= 0 (the same as
    1/2 + atan(z) / pi, without its cancellation in the lower tail; for z >= 0 the
    result is at least 1/2, so the subtraction does not cancel)."""
    with mpmath.extraprec(60):
        z = (x - mu) / b
        return mpmath.acot(-z) / mpmath.pi if z < 0 else 1 - mpmath.acot(z) / mpmath.pi


def _mp_cauchy_sf(mu: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Q(x) = F(mu - (x - mu)): the distribution is symmetric about mu."""
    with mpmath.extraprec(60):
        return _mp_cauchy_cdf(mu, b, 2 * mu - x)


def _cauchy_ppf_grid() -> list[tuple[float, float, float]]:
    # The gamma_ppf probabilities, the doubles next to 1/2 (where the quantile is
    # mu + b tan(pi d) with d tiny) and p = 1/8 and 7/8 (tan = -+(1 + sqrt 2), so
    # with the last pair b tan overflows although the quantile at 7/8 is finite).
    ps = sorted({p for _, p in _gamma_ppf_grid()} | {0.5 - 2.0**-54, 0.5 + 2.0**-53, 0.125, 0.875})
    return [(mu, b, p) for mu, b in _CAUCHY_PARAMETERS for p in ps]


def _cauchy_ppf_tags(mu: float, b: float, p: float) -> list[str]:
    # The quantile is -inf at p = 0 and +inf at p = 1.
    return [POLE] if p in (0.0, 1.0) else []


def _mp_cauchy_ppf(mu: mpmath.mpf, b: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = mu + b tan(pi (p - 1/2)) = mu - b cospi(p) / sinpi(p): no subtraction from p
    (mpmath reduces the arguments of sinpi and cospi exactly, so both are accurate
    next to 0 and 1, and cospi(1/2) is exactly 0)."""
    with mpmath.extraprec(60):
        return mu - b * mpmath.cospi(p) / mpmath.sinpi(p)


def cauchy_ppf_cancellation(mu: float, b: float, p: float) -> mpmath.mpf:
    """(|mu| + |b T|) / |mu + b T| with b T the quantile's tangent part: the ulps
    that adding mu to a rounded b T can cost."""
    m = mpmath.mpf(mu)
    x = _mp_cauchy_ppf(m, mpmath.mpf(b), mpmath.mpf(p))
    return (abs(m) + abs(x - m)) / abs(x)


_CAUCHY_PARAMETERS_DESCRIPTION = (
    "(mu, b) in {(0, 1), (0, 2^-500), (0, 2^500), (-10, 1/4), (100, 8), (-2^1023, 2^1023)}"
)
_CAUCHY_X_DESCRIPTION = (
    f"{_CAUCHY_PARAMETERS_DESCRIPTION} x x = mu + b z for z = -40..40 step 1/4, "
    "z = +-{50, 100, 1e3, 1e4, 1e5, 1e8, 1e10, 1e15, 1e16, 1e17, 1e20, 1e50, 1e100, 1e300, "
    "1e305, 1e308} and z = +-2^-30, +-2^-60, kept where x is an exact, finite double that is "
    "0 or normal (so z = (x - mu) / b is exact)"
)

CAUCHY_CDF = Spec(
    name="cauchy_cdf",
    description=(
        "Cauchy distribution cumulative distribution function 1/2 + atan(z) / pi, "
        "z = (x - mu) / b (location mu, scale b)"
    ),
    inputs=("mu", "b", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, mu and b are exact doubles, b is a power of two and x - mu = b z with z a double, "
        "so z is exact; F's condition number in z, |z F'(z) / F(z)|, is at most 1, so a careful "
        "method (acot(-z) / pi for z < 0, 1 - acot(z) / pi for z >= 0) is within a few ulps. "
        "1e-13 as for the other functions; results below 2.2e-308 are judged relative to that "
        "number (conventions)"
    ),
    grid_description=(
        f"{_CAUCHY_X_DESCRIPTION}; F is about 1 / (pi |z|) in the lower tail (subnormal from "
        "about z = -1.4e307) and rounds to 1 above about z = 5.7e15"
    ),
    grid=_cauchy_x_grid,
    mpmath_fn=_mp_cauchy_cdf,
    scipy_name="scipy.stats.cauchy.cdf",
    scipy_fn=lambda mu, b, x: st.cauchy.cdf(x, loc=mu, scale=b),
)

CAUCHY_SF = Spec(
    name="cauchy_sf",
    description=(
        "Cauchy distribution survival function 1/2 - atan(z) / pi, z = (x - mu) / b "
        "(location mu, scale b)"
    ),
    inputs=("mu", "b", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "the mirror image of cauchy_cdf (Q(x) = F(2 mu - x)), with the same exact z and the "
        "same conditioning. 1e-13 as for the other functions; results below 2.2e-308 are "
        "judged relative to that number (conventions)"
    ),
    grid_description=(
        f"the cauchy_cdf grid: {_CAUCHY_X_DESCRIPTION}; Q rounds to 1 below about "
        "z = -5.7e15 and is about 1 / (pi z) in the upper tail (subnormal from about "
        "z = 1.4e307)"
    ),
    grid=_cauchy_x_grid,
    mpmath_fn=_mp_cauchy_sf,
    scipy_name="scipy.stats.cauchy.sf",
    scipy_fn=lambda mu, b, x: st.cauchy.sf(x, loc=mu, scale=b),
)

CAUCHY_PPF = Spec(
    name="cauchy_ppf",
    description="Cauchy distribution quantile function mu + b tan(pi (p - 1/2))",
    inputs=("mu", "b", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double; a careful method forms T = tan(pi (p - 1/2)) within a few ulps "
        "(tan(pi d) with d = p - 1/2 exact (Sterbenz) for 1/4 <= p <= 3/4, -cot(pi p) below "
        "and cot(pi (1 - p)) above, 1 - p exact), and exactly -1, 0 and 1 at p = 1/4, 1/2 and "
        "3/4; then mu + b T costs (|mu| + |b T|) / |mu + b T| ulps, at most about 8.8 on this "
        "grid where the quantile is not 0 (mu = -10, b = 1/4, p = 0.99; tested), 1e-15, and "
        "nothing where it is 0 (mu = -b, p = 3/4: T = 1 exactly). 1e-13 as for the other "
        "quantiles"
    ),
    grid_description=(
        f"{_CAUCHY_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities (p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1), the doubles next to 1/2 (1/2 - 2^-54, 1/2 + 2^-53) and "
        "p = 1/8, 7/8; p = 0 and 1 are tagged pole (-inf and +inf), quantiles beyond the "
        "largest double are tagged overflow"
    ),
    grid=_cauchy_ppf_grid,
    mpmath_fn=_mp_cauchy_ppf,
    scipy_name="scipy.stats.cauchy.ppf",
    scipy_fn=lambda mu, b, p: st.cauchy.ppf(p, loc=mu, scale=b),
    input_tags=_cauchy_ppf_tags,
)

# --- Levy distribution (location mu, scale c) -----------------------------------

# F(x) = erfc(sqrt(c / (2 (x - mu)))) for x > mu and 0 otherwise. The Gumbel
# pairs, for the same reasons: scales are powers of two and every grid x is
# mu + c t exactly, so t = (x - mu) / c is exact; the last pair puts x - mu
# above the largest double although t is only 2 to 2.75.
_LEVY_PARAMETERS = _GUMBEL_PARAMETERS
# t: below and at the start of the support (F = 0, Q = 1); next to it, where F is
# erfc(s) with s = 1 / sqrt(2 t) large (F is subnormal from about t = 7.1e-4 and 0
# below about t = 6.7e-4); the bulk in steps of 1/4096, 1/64 and 1/4; the heavy
# upper tail, where Q is about sqrt(2 / (pi t)), out to t = 2^1523, which only the
# scale 2^-500 reaches with a finite x (from t = 2^1022, 1 / (2 t) is below the
# smallest normal double and, unless t is a power of two, loses bits; from
# t = 2^1075 it is below every double).
_LEVY_T = tuple(
    sorted(
        {mpmath.mpf(t) for t in (-1.0, -(2.0**-30), 0.0, 2.0**-60, 2.0**-30, 1e-10, 1e-5)}
        | {mpmath.mpf(j) / 65536 for j in range(40, 64)}
        | {mpmath.mpf(j) / 4096 for j in range(4, 64)}
        | {mpmath.mpf(j) / 64 for j in range(1, 64)}
        | {mpmath.mpf(j) / 4 for j in range(4, 161)}
        | {mpmath.mpf(t) for t in (50.0, 100.0, 1e3, 1e4, 1e5, 1e8, 1e10, 1e15, 1e16, 1e17,
                                   1e20, 1e50, 1e100, 1e300, 1e305, 1e308)}
        | {3 * mpmath.mpf(2) ** k for k in (1030, 1050, 1070)}
        | {mpmath.mpf(2) ** k for k in (1100, 1200, 1500, 1523)}
    )
)  # fmt: skip


def _levy_x_grid() -> list[tuple[float, float, float]]:
    # x = mu + c t, kept where it is an exact, finite double that is 0 or normal.
    points = []
    for mu, c in _LEVY_PARAMETERS:
        for t in _LEVY_T:
            # 2200 bits hold mu + c t exactly for every value on this grid.
            with mpmath.workprec(2200):
                exact = mpmath.mpf(mu) + mpmath.mpf(c) * t
                if abs(exact) >= mpmath.mpf(2) ** 1024:
                    continue
                x = float(exact)
                if mpmath.mpf(x) == exact and (x == 0 or abs(x) >= 2.0**-1022):
                    points.append((mu, c, x))
    return sorted(points)


def _levy_s(mu: mpmath.mpf, c: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    # s = sqrt(c / (2 (x - mu))) for x > mu, with 60 extra bits.
    with mpmath.extraprec(60):
        return mpmath.sqrt(c / (2 * (x - mu)))


def _mp_levy_cdf(mu: mpmath.mpf, c: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x) = erfc(s) with s = sqrt(c / (2 (x - mu))) for x > mu, and 0 for x <= mu
    (_mp_erfc: mpmath's erfc beyond s = 1e6 overflows internally)."""
    if x <= mu:
        return mpmath.mpf(0)
    with mpmath.extraprec(60):
        return _mp_erfc(_levy_s(mu, c, x))


def _mp_levy_sf(mu: mpmath.mpf, c: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Q(x) = erf(s) for x > mu (mpmath's erf keeps full relative precision for tiny
    s), and 1 for x <= mu."""
    if x <= mu:
        return mpmath.mpf(1)
    with mpmath.extraprec(60):
        return mpmath.erf(_levy_s(mu, c, x))


def levy_cdf_condition_in_u(mu: float, c: float, x: float) -> mpmath.mpf:
    """|u dF/du / F| at u = c / (2 (x - mu)): F = erfc(sqrt u), so this is
    sqrt(u) e^-u / (sqrt(pi) erfc(sqrt u)), about u + 1/2 for large u."""
    s = _levy_s(*(mpmath.mpf(v) for v in (mu, c, x)))
    return s * mpmath.exp(-(s**2)) / (mpmath.sqrt(mpmath.pi) * _mp_erfc(s))


def levy_sf_condition_in_u(mu: float, c: float, x: float) -> mpmath.mpf:
    """|u dQ/du / Q| at u = c / (2 (x - mu)): Q = erf(sqrt u), so this is
    sqrt(u) e^-u / (sqrt(pi) erf(sqrt u)), at most 1/2 (its limit as u -> 0)."""
    s = _levy_s(*(mpmath.mpf(v) for v in (mu, c, x)))
    return s * mpmath.exp(-(s**2)) / (mpmath.sqrt(mpmath.pi) * mpmath.erf(s))


def _levy_ppf_grid() -> list[tuple[float, float, float]]:
    # The gamma_ppf probabilities and the doubles next to 1/2 and 1.
    ps = sorted({p for _, p in _gamma_ppf_grid()} | {0.5 - 2.0**-54, 0.5 + 2.0**-53, 1 - 2.0**-53})
    return [(mu, c, p) for mu, c in _LEVY_PARAMETERS for p in ps]


def _levy_ppf_tags(mu: float, c: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1. At p = 0 it is mu, the start of the support:
    # a finite value, so not a pole.
    return [POLE] if p == 1.0 else []


def _mp_levy_ppf(mu: mpmath.mpf, c: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """x = mu + c / (2 erfc^-1(p)^2) = mu + c / z^2 with z = Phi^-1(p / 2), as
    erfc^-1(p) = -z / sqrt(2) (p / 2 is exact in mpmath); x = mu at p = 0."""
    if p == 0:
        return mpmath.mpf(mu)
    with mpmath.extraprec(60):
        return mu + c / _mp_norm_ppf(p / 2) ** 2


def levy_ppf_cancellation(mu: float, c: float, p: float) -> mpmath.mpf:
    """(|mu| + |c / z^2|) / |mu + c / z^2|: the ulps that adding mu to a rounded
    c / z^2 can cost."""
    m = mpmath.mpf(mu)
    x = _mp_levy_ppf(m, mpmath.mpf(c), mpmath.mpf(p))
    return (abs(m) + abs(x - m)) / abs(x)


_LEVY_PARAMETERS_DESCRIPTION = (
    "(mu, c) in {(0, 1), (0, 2^-500), (0, 2^500), (-10, 1/4), (100, 8), (-2^1023, 2^1023)}"
)
_LEVY_X_DESCRIPTION = (
    f"{_LEVY_PARAMETERS_DESCRIPTION} x x = mu + c t for t in {{-1, -2^-30, 0, 2^-60, 2^-30, "
    "1e-10, 1e-5}, t = j / 65536 for j = 40..63, j / 4096 for j = 4..63, j / 64 for "
    "j = 1..63, j / 4 for j = 4..160, t in {50, 100, 1e3, 1e4, 1e5, 1e8, 1e10, 1e15, 1e16, "
    "1e17, 1e20, 1e50, 1e100, 1e300, 1e305, 1e308}, t = 3 2^1030, 3 2^1050, 3 2^1070 and "
    "t = 2^1100, 2^1200, 2^1500, 2^1523, "
    "kept where x is an exact, finite double that is 0 or normal (so t = (x - mu) / c is "
    "exact)"
)

LEVY_CDF = Spec(
    name="levy_cdf",
    description=(
        "Levy distribution cumulative distribution function erfc(sqrt(c / (2 (x - mu)))) "
        "(location mu, scale c)"
    ),
    inputs=("mu", "c", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "x, mu and c are exact doubles, c is a power of two and x - mu = c t with t exact on "
        "this grid. F = erfc(sqrt u) with u = c / (2 (x - mu)) has condition number about "
        "u + 1/2 in u, at most 698 at a normal result on this grid (tested), so rounding u "
        "once costs up to 7.8e-14 there (as for gumbel_cdf); a careful method then forms "
        "exp(-u) erfcx(sqrt u), which costs a few ulps more (the square root inside erfc "
        "would double the loss). 1e-13 as for the other functions; results below 2.2e-308 "
        "are judged relative to that number (conventions)"
    ),
    grid_description=(
        f"{_LEVY_X_DESCRIPTION}; F is 0 for t <= 0, subnormal from about t = 7.1e-4, 0 again "
        "(underflow) below about t = 6.7e-4, and about 1 - sqrt(2 / (pi t)) in the upper tail"
    ),
    grid=_levy_x_grid,
    mpmath_fn=_mp_levy_cdf,
    scipy_name="scipy.stats.levy.cdf",
    scipy_fn=lambda mu, c, x: st.levy.cdf(x, loc=mu, scale=c),
)

LEVY_SF = Spec(
    name="levy_sf",
    description=(
        "Levy distribution survival function erf(sqrt(c / (2 (x - mu)))) (location mu, scale c)"
    ),
    inputs=("mu", "c", "x"),
    tolerance=1e-13,
    tolerance_reason=(
        "t = (x - mu) / c is exact on this grid (see levy_cdf). Q = erf(sqrt u) with "
        "u = c / (2 (x - mu)) has condition number at most 1/2 in u, so rounding u and the erf "
        "cost a few ulps; where u is below the smallest normal double (t = 1e308, and "
        "t >= 3 2^1030 with scale 2^-500) a careful method forms sqrt(c / 2) / sqrt(x - mu) "
        "instead. 1e-13 as "
        "for the other functions; results below 2.2e-308 are judged relative to that number "
        "(conventions)"
    ),
    grid_description=(
        f"the levy_cdf grid: {_LEVY_X_DESCRIPTION}; Q is 1 for t <= 0, rounds to 1 below about "
        "t = 0.014 and is about sqrt(2 / (pi t)) in the upper tail (down to 4.7e-230 at "
        "t = 2^1523)"
    ),
    grid=_levy_x_grid,
    mpmath_fn=_mp_levy_sf,
    scipy_name="scipy.stats.levy.sf",
    scipy_fn=lambda mu, c, x: st.levy.sf(x, loc=mu, scale=c),
)

LEVY_PPF = Spec(
    name="levy_ppf",
    description="Levy distribution quantile function mu + c / (2 erfc^-1(p)^2)",
    inputs=("mu", "c", "p"),
    tolerance=1e-13,
    tolerance_reason=(
        "p is an exact double; a careful method forms erfc^-1(p) within a few ulps (Boost's "
        "erf_inv reaches about 1e-16; 1 - p is exact for p >= 1/2), then c / (2 y^2) costs a "
        "few ulps more and mu + c / (2 y^2) costs (|mu| + |c / (2 y^2)|) / |x| ulps, at most "
        "about 7.2 on this grid (tested), 8e-16; it is exactly mu at p = 0. 1e-13 as for the "
        "other quantiles"
    ),
    grid_description=(
        f"{_LEVY_PARAMETERS_DESCRIPTION} x the gamma_ppf probabilities (p = 0, 5e-324, "
        "1e-300, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, "
        "1/2 - 1e-10, 1/2, 1/2 + 1e-10, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15 and 1) and the doubles next to 1/2 (1/2 - 2^-54, 1/2 + 2^-53) "
        "and 1 (1 - 2^-53); the quantile is mu at p = 0, p = 1 is tagged pole (+inf), "
        "quantiles beyond the largest double are tagged overflow"
    ),
    grid=_levy_ppf_grid,
    mpmath_fn=_mp_levy_ppf,
    scipy_name="scipy.stats.levy.ppf",
    scipy_fn=lambda mu, c, p: st.levy.ppf(p, loc=mu, scale=c),
    input_tags=_levy_ppf_tags,
)

# --- hypothesis tests: one-sample t-test --------------------------------------

# Deterministic samples x_i = loc + z_i with z_i = ((7919 i) mod 1000 - 500) / 1000
# (a fixed scrambled pattern of exact decimals, sd about 0.29; no randomness),
# and the hypothesised mean mu = loc + d, so the statistic t is about
# -d sqrt(n) / 0.29: from 0.003 to 3400 in size.
_TTEST_N = (2, 3, 5, 10, 30, 100)
_TTEST_LOC = (0.0, 10.0)
_TTEST_D = (-30.0, -3.0, -0.3, -0.01, 0.001, 0.1, 1.0, 10.0, 100.0)
# The alternative hypothesis, coded as a number: -1 less, 0 two-sided, 1 greater.
_TTEST_ALTERNATIVES = {-1.0: "less", 0.0: "two-sided", 1.0: "greater"}


def _ttest_sample(n: int, loc: float) -> tuple[float, ...]:
    return tuple(loc + ((7919 * i) % 1000 - 500) / 1000 for i in range(n))


def _mp_ttest_parts(x: list[mpmath.mpf], mu: mpmath.mpf):
    """The exact mean, variance (n - 1 denominator) and statistic t."""
    n = len(x)
    with mpmath.extraprec(60):
        mean = mpmath.fsum(x) / n
        var = mpmath.fsum((v - mean) ** 2 for v in x) / (n - 1)
        t = (mean - mu) / mpmath.sqrt(var / n)
    return mean, var, t


def ttest_condition(x: tuple[float, ...], mu: float, alternative: float) -> tuple:
    """(kappa_t, kappa_p): the relative condition numbers of the statistic and of
    the p-value as a function of t.

    kappa_t bounds |dt / t| per relative perturbation eps of every x_i and mu:
    (mean |x_i| + |mu|) / |mean - mu| from the numerator and
    sum |x_i - mean| |x_i| / ((n - 1) s^2) from the standard deviation.
    kappa_p = |t f(t) / P| with P the p-value's tail (both tails at |t| for the
    two-sided test, whose p-value is 2 sf(|t|)).
    """
    xs = [mpmath.mpf(v) for v in x]
    mu_m = mpmath.mpf(mu)
    n = len(xs)
    mean, var, t = _mp_ttest_parts(xs, mu_m)
    kappa_t = (mpmath.fsum(abs(v) for v in xs) / n + abs(mu_m)) / abs(mean - mu_m)
    kappa_t += mpmath.fsum(abs(v - mean) * abs(v) for v in xs) / ((n - 1) * var)
    df = mpmath.mpf(n - 1)
    if alternative == 0:
        tail = _mp_t_sf(df, abs(t))
    elif alternative > 0:
        tail = _mp_t_sf(df, t)
    else:
        tail = _mp_t_cdf(df, t)
    kappa_p = abs(t * _mp_t_pdf(df, t) / tail)
    return kappa_t, kappa_p


def ttest_rounding_bound(x: tuple[float, ...], mu: float, alternative: float) -> float:
    """What a careful double-precision method may lose in the p-value: the statistic
    from sums of n rounded terms is off by up to about (kappa_t + n) 2^-53
    relative, and the p-value by kappa_p times that."""
    kappa_t, kappa_p = ttest_condition(x, mu, alternative)
    return float(max(kappa_p, 1) * (kappa_t + len(x)) * mpmath.mpf(2) ** -53)


def _ttest_grid() -> list[tuple[tuple[float, ...], float, float]]:
    points = []
    for n in _TTEST_N:
        for loc in _TTEST_LOC:
            x = _ttest_sample(n, loc)
            for d in _TTEST_D:
                points.extend((x, loc + d, alt) for alt in _TTEST_ALTERNATIVES)
    return points


def _mp_ttest_1samp(x: list[mpmath.mpf], mu: mpmath.mpf, alternative: mpmath.mpf):
    """(t, p) for the one-sample t-test with n - 1 degrees of freedom."""
    _, _, t = _mp_ttest_parts(x, mu)
    df = mpmath.mpf(len(x) - 1)
    if alternative == 0:
        p = 2 * _mp_t_sf(df, abs(t))
    elif alternative > 0:
        p = _mp_t_sf(df, t)
    else:
        p = _mp_t_cdf(df, t)
    return t, p


def _sc_ttest_1samp(x: list[float], mu: float, alternative: float) -> tuple[float, float]:
    r = st.ttest_1samp(x, mu, alternative=_TTEST_ALTERNATIVES[alternative])
    return float(r.statistic), float(r.pvalue)


TTEST_1SAMP = Spec(
    name="ttest_1samp",
    description=(
        "One-sample t-test: statistic t = (mean - mu) / (s / sqrt(n)) and its p-value from "
        "Student's t with n - 1 degrees of freedom; alternative -1 less, 0 two-sided, 1 greater"
    ),
    inputs=("x", "mu", "alternative"),
    outputs=("statistic", "pvalue"),
    tolerance=5e-12,
    tolerance_reason=(
        "For both outputs. The inputs are exact, but t is formed from sums of n rounded terms, "
        "so a careful method is off by up to about (kappa_t + n) 2^-53, kappa_t being t's "
        "relative condition number in the data and mu (ttest_condition); the p-value adds the "
        "factor kappa_p = |t f(t) / P| (at most 99 here). The largest bound "
        "max(kappa_p, 1) (kappa_t + n) 2^-53 on this grid is 4.5e-12 (n = 100, loc = 10, "
        "mu = 9.99, where t = 0.017 and kappa_t = 4.0e4: mean - mu cancels); 5e-12 is the "
        "smallest round value above it (tested). SciPy's worst error is 1.7e-12, at that point"
    ),
    grid_description=(
        "samples x_i = loc + ((7919 i) mod 1000 - 500) / 1000, i = 0..n-1, for n in {2, 3, 5, "
        "10, 30, 100} and loc in {0, 10}; mu = loc + d for d in {-30, -3, -0.3, -0.01, 0.001, "
        "0.1, 1, 10, 100}; each with alternative -1 (less), 0 (two-sided) and 1 (greater). |t| "
        "runs from 0.017 to 3410 and the p-value from 4.3e-253 to 1"
    ),
    grid=_ttest_grid,
    mpmath_fn=_mp_ttest_1samp,
    scipy_name="scipy.stats.ttest_1samp",
    scipy_fn=_sc_ttest_1samp,
)

# --- hypothesis tests: chi-squared goodness of fit -------------------------------

# Expected frequencies e_i either equal (total / k, statrs's default) or
# proportional to 1, 2, .., k, and observed counts o_i = round(e_i (1 + delta z_i))
# with the z pattern of the t-test samples; the expected frequencies are then
# rescaled to the counts' sum, as the test requires. The statistic is about
# delta^2 total / 12: from about 1e-5 to 900.
_CHISQ_K = (2, 3, 5, 10, 30)
_CHISQ_TOTAL = (100, 10**4, 10**6)
_CHISQ_DELTA = (0.001, 0.01, 0.03, 0.1)


def _chisquare_grid() -> list[tuple[tuple[float, ...], tuple[float, ...], float]]:
    points = []
    for k in _CHISQ_K:
        z = [((7919 * i) % 1000 - 500) / 1000 for i in range(k)]
        for weights in ([1] * k, list(range(1, k + 1))):
            for total in _CHISQ_TOTAL:
                for delta in _CHISQ_DELTA:
                    e0 = [total * w / sum(weights) for w in weights]
                    obs = tuple(float(round(e * (1 + delta * zi))) for e, zi in zip(e0, z))
                    exp = tuple(sum(obs) * w / sum(weights) for w in weights)
                    if obs == exp:
                        continue  # rounding left the counts at their expectation
                    for ddof in (0.0, 1.0):
                        # A small total can round two deltas to the same counts.
                        if ddof < k - 1 and (obs, exp, ddof) not in points:
                            points.append((obs, exp, ddof))
    return points


def _mp_chisquare(obs: list[mpmath.mpf], exp: list[mpmath.mpf], ddof: mpmath.mpf):
    """(statistic, p-value): sum (o - e)^2 / e and Q(dof / 2, statistic / 2) with
    dof = k - 1 - ddof (the chi-squared survival function, the gamma_sf route)."""
    with mpmath.extraprec(60):
        stat = mpmath.fsum((o - e) ** 2 / e for o, e in zip(obs, exp, strict=True))
        dof = len(obs) - 1 - ddof
        return stat, _mp_gamma_sf(dof / 2, stat / 2)


def chisquare_rounding_bound(obs: tuple, exp: tuple, ddof: float) -> float:
    """max(kappa_p, 1) (kappa_s + k) 2^-53: what a careful method may lose.

    kappa_s = sum |o - e| |o + e| / e / statistic is the statistic's relative
    condition number in the expected frequencies (the counts are exact);
    kappa_p = |x f(x) / Q| is the p-value's in the statistic x."""
    o = [mpmath.mpf(v) for v in obs]
    e = [mpmath.mpf(v) for v in exp]
    stat, _ = _mp_chisquare(o, e, mpmath.mpf(ddof))
    kappa_s = mpmath.fsum(abs(a - b) * abs(a + b) / b for a, b in zip(o, e, strict=True)) / stat
    a = (len(o) - 1 - ddof) / mpmath.mpf(2)
    kappa_p = gamma_sf_condition(a, stat / 2)
    return float(max(kappa_p, 1) * (kappa_s + len(o)) * mpmath.mpf(2) ** -53)


def _sc_chisquare(obs: list[float], exp: list[float], ddof: float) -> tuple[float, float]:
    r = st.chisquare(obs, exp, ddof=int(ddof))
    return float(r.statistic), float(r.pvalue)


CHISQUARE = Spec(
    name="chisquare",
    description=(
        "Chi-squared goodness-of-fit test: statistic sum (o - e)^2 / e over k categories and "
        "its p-value, the chi-squared survival function with k - 1 - ddof degrees of freedom"
    ),
    inputs=("f_obs", "f_exp", "ddof"),
    outputs=("statistic", "pvalue"),
    tolerance=5e-12,
    tolerance_reason=(
        "For both outputs, as for ttest_1samp. The counts are exact and the statistic is a "
        "sum of k positive terms, off by up to about (kappa_s + k) 2^-53 in a careful method, "
        "kappa_s being its relative condition number in the expected frequencies "
        "(chisquare_rounding_bound); the p-value adds kappa_p = |x f(x) / Q|, about half the "
        "statistic in the tail (431 at most where the p-value is a normal double). The largest "
        "bound max(kappa_p, 1) (kappa_s + k) 2^-53 there is 4.2e-12 (k = 30 equal expected "
        "frequencies, total 1e6, delta 0.1, ddof 1: statistic 888, p 6.1e-169); 5e-12 as for "
        "ttest_1samp (tested). P-values that round to 0 (4 points) are exact in any method. "
        "SciPy's worst error is 1.9e-13"
    ),
    grid_description=(
        "k in {2, 3, 5, 10, 30} categories; expected frequencies equal or proportional to 1, 2, "
        "..., k; totals 100, 1e4, 1e6; observed counts round(e_i (1 + delta z_i)) with "
        "z_i = ((7919 i) mod 1000 - 500) / 1000 and delta in {0.001, 0.01, 0.03, 0.1} (skipped "
        "where every count equals its expectation), the expected frequencies rescaled to the "
        "counts' sum; ddof 0 and, where k > 2, 1. Statistic from 1e-3 to 2120, p-value from 0 "
        "to 1"
    ),
    grid=_chisquare_grid,
    mpmath_fn=_mp_chisquare,
    scipy_name="scipy.stats.chisquare",
    scipy_fn=_sc_chisquare,
)

# --- hypothesis tests: one-way analysis of variance -----------------------------

# Groups of the given sizes, all values in one sample x with the sizes as a
# second input; value i (counted over all groups) of group j is
# loc + z_i + delta j, with the t-test pattern z.
_ANOVA_SIZES = ((2, 2), (3, 3, 3), (2, 4, 6, 8, 10), (30, 30, 30))
_ANOVA_LOC = (0.0, 100.0)
_ANOVA_DELTA = (0.001, 0.03, 0.3, 1.0, 3.0, 30.0)


def _f_oneway_grid() -> list[tuple[tuple[float, ...], tuple[float, ...]]]:
    points = []
    for sizes in _ANOVA_SIZES:
        for loc in _ANOVA_LOC:
            for delta in _ANOVA_DELTA:
                x, i = [], 0
                for j, size in enumerate(sizes):
                    for _ in range(size):
                        x.append(loc + ((7919 * i) % 1000 - 500) / 1000 + delta * j)
                        i += 1
                points.append((tuple(x), tuple(float(v) for v in sizes)))
    return points


def _anova_groups(x: list, sizes: list) -> list[list]:
    groups, start = [], 0
    for size in sizes:
        groups.append(x[start : start + int(size)])
        start += int(size)
    return groups


def _mp_anova_parts(x: list[mpmath.mpf], sizes: list[mpmath.mpf]):
    """(ssb, ssw, d1, d2) with 60 extra bits."""
    groups = _anova_groups(x, sizes)
    with mpmath.extraprec(60):
        mean = mpmath.fsum(x) / len(x)
        means = [mpmath.fsum(g) / len(g) for g in groups]
        ssb = mpmath.fsum(len(g) * (m - mean) ** 2 for g, m in zip(groups, means, strict=True))
        ssw = mpmath.fsum((v - m) ** 2 for g, m in zip(groups, means, strict=True) for v in g)
    return ssb, ssw, mpmath.mpf(len(groups) - 1), mpmath.mpf(len(x) - len(groups))


def _mp_f_oneway(x: list[mpmath.mpf], sizes: list[mpmath.mpf]):
    """(F, p-value): F = (ssb / (k - 1)) / (ssw / (N - k)) and the F survival function."""
    ssb, ssw, d1, d2 = _mp_anova_parts(x, sizes)
    with mpmath.extraprec(60):
        f = (ssb / d1) / (ssw / d2)
    return f, _mp_f_sf(d1, d2, f)


def f_oneway_rounding_bound(x: tuple, sizes: tuple) -> float:
    """max(kappa_p, 1) (kappa_F + N) 2^-53: what a careful (two-pass) method may lose.

    kappa_F = kappa_ssb + kappa_ssw bounds F's relative condition number in the
    data: perturbing x_i by eps x_i moves ssw by up to 2 eps |x_i - m_j| |x_i|
    and ssb by up to 2 eps |m_j - m| |x_i| (1 - n_j / N); kappa_p = |F f(F) / Q|.
    """
    xs = [mpmath.mpf(v) for v in x]
    groups = _anova_groups(xs, sizes)
    n = len(xs)
    ssb, ssw, d1, d2 = _mp_anova_parts(xs, [mpmath.mpf(v) for v in sizes])
    mean = mpmath.fsum(xs) / n
    k_ssw = k_ssb = mpmath.mpf(0)
    for g in groups:
        m = mpmath.fsum(g) / len(g)
        k_ssw += mpmath.fsum(2 * abs(v - m) * abs(v) for v in g)
        k_ssb += mpmath.fsum(2 * abs(m - mean) * abs(v) * (1 - mpmath.mpf(len(g)) / n) for v in g)
    kappa_f = k_ssw / ssw + k_ssb / ssb
    f = (ssb / d1) / (ssw / d2)
    kappa_p = f_condition(float(d1), float(d2), float(f), upper=True)
    return float(max(kappa_p, 1) * (kappa_f + n) * mpmath.mpf(2) ** -53)


def _sc_f_oneway(x: list[float], sizes: list[float]) -> tuple[float, float]:
    r = st.f_oneway(*_anova_groups(x, sizes))
    return float(r.statistic), float(r.pvalue)


F_ONEWAY = Spec(
    name="f_oneway",
    description=(
        "One-way analysis of variance: F = (between-group mean square) / (within-group mean "
        "square) for k groups of N values and its p-value, the F survival function with "
        "k - 1 and N - k degrees of freedom; the groups are given as one sample x and their "
        "sizes"
    ),
    inputs=("x", "sizes"),
    outputs=("statistic", "pvalue"),
    tolerance=5e-12,
    tolerance_reason=(
        "For both outputs, as for ttest_1samp. The inputs are exact, but F is formed from sums "
        "of N rounded terms, so a careful (two-pass) method is off by up to about "
        "(kappa_F + N) 2^-53, kappa_F being F's relative condition number in the data "
        "(f_oneway_rounding_bound); the p-value adds kappa_p = |F f(F) / Q|. The largest bound "
        "max(kappa_p, 1) (kappa_F + N) 2^-53 on this grid is 4.2e-12 (three groups of 30 "
        "around 100 with delta 30: F = 3.0e5, kappa_p = 43); 5e-12 as for ttest_1samp "
        "(tested). SciPy's worst error at the points it gets within tolerance is 5.0e-12; it "
        "is off at 4: one statistic (5.1e-12, F = 1.2e5) and four p-values below 1e-52 "
        "(7.4e-12 to 7.5e-11)"
    ),
    grid_description=(
        "group sizes (2, 2), (3, 3, 3), (2, 4, 6, 8, 10), (30, 30, 30); value i (counted over "
        "all groups) of group j is loc + ((7919 i) mod 1000 - 500) / 1000 + delta j for loc in "
        "{0, 100} and delta in {0.001, 0.03, 0.3, 1, 3, 30}; x lists all values and sizes the "
        "group sizes. F from 0.33 to 3.0e5, p-value from 5.8e-168 to 0.72"
    ),
    grid=_f_oneway_grid,
    mpmath_fn=_mp_f_oneway,
    scipy_name="scipy.stats.f_oneway",
    scipy_fn=_sc_f_oneway,
)

# --- hypothesis tests: D'Agostino's skewness test ----------------------------

# Samples x_i = loc + u_i + s u_i^2 (+ an outlier) with u_i the t-test pattern
# ((7919 i) mod 1000 - 500) / 1000: s skews them; an outlier of 50 at the first
# value skews them strongly. Plus exactly symmetric samples, whose sample
# skewness is 0.
_SKEW_N = (8, 10, 20, 50, 100)
_SKEW_S = (-3.0, -0.3, 0.01, 0.3, 3.0)


def _skewtest_samples() -> list[tuple[float, ...]]:
    samples = []
    for n in _SKEW_N:
        u = [((7919 * i) % 1000 - 500) / 1000 for i in range(n)]
        for loc in (0.0, 10.0):
            samples.extend(tuple(loc + v + s * v * v for v in u) for s in _SKEW_S)
            samples.append(tuple(loc + v + (50.0 if i == 0 else 0.0) for i, v in enumerate(u)))
            # Symmetric: the values 1..n shifted to loc (their mean, minus (n + 1) / 2).
            samples.append(tuple(loc + i - (n + 1) / 2 for i in range(1, n + 1)))
    return samples


def _skewtest_grid() -> list[tuple[tuple[float, ...], float]]:
    return [(x, alt) for x in _skewtest_samples() for alt in _TTEST_ALTERNATIVES]


def _mp_skew_root_b1(x: list[mpmath.mpf]) -> mpmath.mpf:
    """The sample skewness m3 / m2^(3/2) with the biased (1/n) central moments."""
    n = len(x)
    mean = mpmath.fsum(x) / n
    m2 = mpmath.fsum((v - mean) ** 2 for v in x) / n
    m3 = mpmath.fsum((v - mean) ** 3 for v in x) / n
    return m3 / m2 ** mpmath.mpf(1.5)


def _mp_skewtest_z(x: list[mpmath.mpf]) -> mpmath.mpf:
    """D'Agostino's transform of the sample skewness to a standard normal z
    (D'Agostino, Belanger and D'Agostino 1990, eqs. 4-9), as in SciPy:
    y = b1 sqrt((n + 1)(n + 3) / (6 (n - 2))), W^2 = -1 + sqrt(2 (beta2 - 1)),
    delta = 1 / sqrt(log W), alpha = sqrt(2 / (W^2 - 1)), z = delta asinh(y / alpha).
    z = 0 for a sample with zero skewness."""
    n = mpmath.mpf(len(x))
    with mpmath.extraprec(60):
        y = _mp_skew_root_b1(x) * mpmath.sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
        beta2 = (
            3 * (n**2 + 27 * n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
        )
        w2 = -1 + mpmath.sqrt(2 * (beta2 - 1))
        delta = 1 / mpmath.sqrt(mpmath.log(w2) / 2)
        alpha = mpmath.sqrt(2 / (w2 - 1))
        return delta * mpmath.asinh(y / alpha)


def _mp_norm_tail(z: mpmath.mpf, alternative: float) -> mpmath.mpf:
    """The p-value of a standard normal z: 2 Q(|z|), Q(z) or Phi(z)."""
    if alternative == 0:
        return mpmath.erfc(abs(z) / mpmath.sqrt(2))
    if alternative > 0:
        return mpmath.erfc(z / mpmath.sqrt(2)) / 2
    return mpmath.erfc(-z / mpmath.sqrt(2)) / 2


def _mp_skewtest(x: list[mpmath.mpf], alternative: mpmath.mpf):
    z = _mp_skewtest_z(x)
    return z, _mp_norm_tail(z, alternative)


def skewtest_rounding_bound(x: tuple[float, ...], alternative: float) -> float:
    """max(kappa_p, 1) (kappa_z + n) 2^-53: what a careful method may lose.

    kappa_z = kappa_b1 times the asinh step's |u / (sqrt(1 + u^2) asinh u)| <= 1
    (u = y / alpha), with kappa_b1 = sum |x_i d b1 / d x_i| / |b1| from
    d m2 / d x_i = 2 d_i / n and d m3 / d x_i = 3 (d_i^2 - m2) / n (d = x - mean);
    kappa_p = |z phi(z) / P|."""
    xs = [mpmath.mpf(v) for v in x]
    n = len(xs)
    mean = mpmath.fsum(xs) / n
    d = [v - mean for v in xs]
    m2 = mpmath.fsum(v**2 for v in d) / n
    m3 = mpmath.fsum(v**3 for v in d) / n
    if m3 == 0:
        return 0.0  # the symmetric samples: exact in double arithmetic (see the tolerance reason)
    b1 = m3 / m2 ** mpmath.mpf(1.5)
    grads = [(3 * (di**2 - m2) / n) / m2**1.5 - 1.5 * m3 / m2**2.5 * (2 * di / n) for di in d]
    kappa_b1 = mpmath.fsum(abs(v * g) for v, g in zip(xs, grads, strict=True)) / abs(b1)
    z = _mp_skewtest_z(xs)
    nn = mpmath.mpf(n)
    y = b1 * mpmath.sqrt((nn + 1) * (nn + 3) / (6 * (nn - 2)))
    beta2 = (
        3
        * (nn**2 + 27 * nn - 70)
        * (nn + 1)
        * (nn + 3)
        / ((nn - 2) * (nn + 5) * (nn + 7) * (nn + 9))
    )
    w2 = -1 + mpmath.sqrt(2 * (beta2 - 1))
    u = y / mpmath.sqrt(2 / (w2 - 1))
    kappa_z = kappa_b1 * abs(u / (mpmath.sqrt(1 + u * u) * mpmath.asinh(u)))
    pdf = mpmath.npdf(z)
    kappa_p = abs(z * pdf / _mp_norm_tail(z, alternative)) * (1 if alternative != 0 else 2)
    return float(max(kappa_p, 1) * (kappa_z + n) * mpmath.mpf(2) ** -53)


def _sc_skewtest(x: list[float], alternative: float) -> tuple[float, float]:
    r = st.skewtest(x, alternative=_TTEST_ALTERNATIVES[alternative])
    return float(r.statistic), float(r.pvalue)


SKEWTEST = Spec(
    name="skewtest",
    description=(
        "D'Agostino's skewness test: z from the sample skewness and its p-value from the "
        "standard normal; alternative -1 less, 0 two-sided, 1 greater"
    ),
    inputs=("x", "alternative"),
    outputs=("statistic", "pvalue"),
    tolerance=2e-12,
    tolerance_reason=(
        "For both outputs, by the rule of ttest_1samp. The skewness b1 = m3 / m2^(3/2) is "
        "formed from sums of n rounded terms, so a careful method is off by up to about "
        "(kappa_z + n) 2^-53, kappa_z being z's relative condition number in the data "
        "(skewtest_rounding_bound); the p-value adds kappa_p = |z phi(z) / P|. The largest "
        "bound max(kappa_p, 1) (kappa_z + n) 2^-53 on this grid is 1.6e-12 (n = 100 with an "
        "outlier, z = 12.1); 2e-12 is the smallest round value above it (tested). For the "
        "exactly symmetric samples every operation in the moments is exact in double "
        "arithmetic, so a careful method returns z = 0 exactly (tested). SciPy's worst error "
        "is 1.5e-12 except at those symmetric samples, where it returns z = 1.01 to 1.04 "
        "instead of 0"
    ),
    grid_description=(
        "samples of n in {8, 10, 20, 50, 100} values: loc + u_i + s u_i^2 with "
        "u_i = ((7919 i) mod 1000 - 500) / 1000, loc in {0, 10} and s in {-3, -0.3, 0.01, 0.3, "
        "3}; loc + u_i plus 50 at i = 0 (an outlier); and the exactly symmetric loc + i - "
        "(n + 1) / 2, i = 1..n; each with alternative -1 (less), 0 (two-sided) and 1 "
        "(greater). |z| from 0 to 12.1, p-value from 4.3e-34 to 1"
    ),
    grid=_skewtest_grid,
    mpmath_fn=_mp_skewtest,
    scipy_name="scipy.stats.skewtest",
    scipy_fn=_sc_skewtest,
)

# --- hypothesis tests: Mann-Whitney U test -------------------------------------

# Two samples from the distinct integers v_i = (7919 i) mod 1000: x is v_0..v_{n1-1}
# shifted by d, y is the next n2. Half-integer shifts keep every value distinct (the
# exact method needs that); the tied samples round v_i down to a multiple of 200 and
# shift by multiples of 200, so values repeat within and across the samples.
# method: 0 automatic, 1 exact, 2 asymptotic with and 3 without the continuity
# correction. The exact sizes keep statrs's enumeration of all C(n1 + n2, min(n1, n2))
# splits below 3.2 million.
_MWU_EXACT_SIZES = ((3, 6), (6, 3), (4, 5), (5, 4), (8, 8), (8, 20), (20, 8))
_MWU_LARGE_SIZES = ((9, 9), (10, 25), (25, 10))
_MWU_TIED_SIZES = ((4, 5), (8, 8), (10, 25))
_MWU_SHIFTS = (-1000.5, -200.5, 0.5, 150.5, 1000.5)
_MWU_TIED_SHIFTS = (-1000.0, -200.0, 0.0, 200.0, 1000.0)
_MWU_METHODS = {0.0: "auto", 1.0: "exact", 2.0: "asymptotic", 3.0: "asymptotic"}


def _mwu_samples(n1: int, n2: int, shift: float, tied: bool):
    v = [(7919 * i) % 1000 for i in range(n1 + n2)]
    if tied:
        v = [200 * (w // 200) for w in v]
    return tuple(float(w) + shift for w in v[:n1]), tuple(float(w) for w in v[n1:])


def _mwu_grid() -> list[tuple[tuple[float, ...], tuple[float, ...], float, float]]:
    grid = []
    for sizes, shifts, tied, methods in (
        (_MWU_EXACT_SIZES, _MWU_SHIFTS, False, (0.0, 1.0, 2.0, 3.0)),
        (_MWU_LARGE_SIZES, _MWU_SHIFTS, False, (0.0, 2.0, 3.0)),
        (_MWU_TIED_SIZES, _MWU_TIED_SHIFTS, True, (0.0, 2.0, 3.0)),
    ):
        for n1, n2 in sizes:
            for d in shifts:
                x, y = _mwu_samples(n1, n2, d, tied)
                grid.extend((x, y, m, alt) for m in methods for alt in _TTEST_ALTERNATIVES)
    return grid


def _mwu_twice_ranks(xy: list) -> list[int]:
    """Twice the average ranks (1-based): integers, so everything stays exact."""
    order = sorted(xy)
    first = {}
    count = {}
    for i, v in enumerate(order):
        first.setdefault(v, i + 1)
        count[v] = count.get(v, 0) + 1
    return [2 * first[v] + count[v] - 1 for v in xy]


def _mwu_twice_u1(x: list, y: list) -> int:
    """2 U1 = 2 R1 - n1 (n1 + 1)."""
    n1 = len(x)
    return sum(_mwu_twice_ranks(list(x) + list(y))[:n1]) - n1 * (n1 + 1)


def _mwu_null_counts(n1: int, n2: int) -> list[int]:
    """c[u] = number of the C(n1 + n2, n1) splits with U = u (no ties): the
    coefficients of the Gaussian binomial [n1 + n2, n1]_q, from the recurrence
    c(m, n; u) = c(m - 1, n; u - n) + c(m, n - 1; u) (Mann and Whitney 1947, eq. 2)."""
    # table[m] holds the counts for (m, n) as n grows; start with n = 0: U = 0 only.
    table = [[1] for _ in range(n1 + 1)]
    for n in range(1, n2 + 1):
        new = [[1]]
        for m in range(1, n1 + 1):
            up = table[m]  # (m, n - 1)
            left = new[m - 1]  # (m - 1, n), shifted by n
            c = [0] * (m * n + 1)
            for u, k in enumerate(up):
                c[u] += k
            for u, k in enumerate(left):
                c[u + n] += k
            new.append(c)
        table = new
    return table[n1]


def _mwu_exact_sf(twice_u: int, n1: int, n2: int) -> tuple[int, int]:
    """P(U >= u) under the null hypothesis as (count, total), exactly (no ties, so U
    is a whole number)."""
    if twice_u % 2:
        raise ValueError(f"exact Mann-Whitney p-value at U = {twice_u}/2: U must be whole (ties?)")
    c = _mwu_null_counts(n1, n2)
    return sum(c[twice_u // 2 :]), sum(c)


def _mwu_tie_term(x: list, y: list) -> int:
    counts = {}
    for v in list(x) + list(y):
        counts[v] = counts.get(v, 0) + 1
    return sum(t**3 - t for t in counts.values())


def _mwu_z(twice_u: int, n1: int, n2: int, tie_term: int, continuity: bool) -> mpmath.mpf:
    """(U - n1 n2 / 2 - 0.5) / s with s^2 = n1 n2 / 12 ((n + 1) - sum(t^3 - t) / (n (n - 1))),
    the variance of U with ties (Lehmann 1975, eq. 1.45); without the 0.5 if not continuity.
    Computed as (2 U - n1 n2 - 1) / (2 s) with 4 s^2 = n1 n2 ((n + 1) n (n - 1) - tie_term)
    / (3 n (n - 1)), integers until the square root."""
    n = n1 + n2
    numerator = twice_u - n1 * n2 - (1 if continuity else 0)
    four_var = mpmath.mpf(n1 * n2 * ((n + 1) * n * (n - 1) - tie_term)) / (3 * n * (n - 1))
    return numerator / mpmath.sqrt(four_var)


def _mwu_tested_u(twice_u1: int, n1: int, n2: int, alternative: float) -> tuple[int, int]:
    """2 U for the U the p-value is taken at, and the factor f: U1 (greater), U2
    (less), max(U1, U2) with f = 2 (two-sided), as in SciPy."""
    twice_u2 = 2 * n1 * n2 - twice_u1
    if alternative > 0:
        return twice_u1, 1
    if alternative < 0:
        return twice_u2, 1
    return max(twice_u1, twice_u2), 2


def _mwu_ties(x: list, y: list) -> bool:
    return len(set(x) | set(y)) < len(x) + len(y)


def _mwu_method(x: list, y: list, method: float) -> float:
    """Method 0 resolved by SciPy's and statrs's rule."""
    if method != 0:
        return method
    return 2.0 if (len(x) > 8 and len(y) > 8) or _mwu_ties(x, y) else 1.0


def _mp_mannwhitneyu(x: list, y: list, method: mpmath.mpf, alternative: mpmath.mpf):
    xs, ys = [float(v) for v in x], [float(v) for v in y]  # the inputs are doubles
    n1, n2 = len(xs), len(ys)
    twice_u1 = _mwu_twice_u1(xs, ys)
    twice_u, f = _mwu_tested_u(twice_u1, n1, n2, float(alternative))
    method = _mwu_method(xs, ys, float(method))
    if method == 1:
        count, total = _mwu_exact_sf(twice_u, n1, n2)
        p = min(mpmath.mpf(f * count) / total, mpmath.mpf(1))
    else:
        z = _mwu_z(twice_u, n1, n2, _mwu_tie_term(xs, ys), continuity=method == 2)
        p = min(f * mpmath.erfc(z / mpmath.sqrt(2)) / 2, mpmath.mpf(1))
    return mpmath.mpf(twice_u1) / 2, p


def mannwhitneyu_rounding_bound(x: tuple, y: tuple, method: float, alternative: float) -> float:
    """max(kappa_p, 1) (kappa_U + n) 2^-53 with kappa_U = 0: U is a sum of
    half-integer ranks, exact in double arithmetic, so U and z move only with the
    few rounded operations in s, counted generously by n; kappa_p = |z phi(z) / Q|
    for the asymptotic methods and 1 for the exact one (a ratio of exact counts)."""
    n1, n2 = len(x), len(y)
    method = _mwu_method(x, y, method)
    kappa_p = mpmath.mpf(1)
    if method != 1:
        twice_u, _ = _mwu_tested_u(_mwu_twice_u1(x, y), n1, n2, alternative)
        z = _mwu_z(twice_u, n1, n2, _mwu_tie_term(x, y), continuity=method == 2)
        kappa_p = max(abs(z * mpmath.npdf(z) / (mpmath.erfc(z / mpmath.sqrt(2)) / 2)), 1)
    return float(kappa_p * (n1 + n2) * mpmath.mpf(2) ** -53)


def _sc_mannwhitneyu(
    x: list[float], y: list[float], method: float, alternative: float
) -> tuple[float, float]:
    r = st.mannwhitneyu(
        x,
        y,
        use_continuity=method != 3,
        alternative=_TTEST_ALTERNATIVES[alternative],
        method=_MWU_METHODS[method],
    )
    return float(r.statistic), float(r.pvalue)


MANNWHITNEYU = Spec(
    name="mannwhitneyu",
    description=(
        "Mann-Whitney U test: U1 = R1 - n1 (n1 + 1) / 2 from the average ranks of x in "
        "x and y, and its p-value; method 0 automatic (exact if n1 <= 8 or n2 <= 8 and no "
        "ties, else 2), 1 exact null distribution, 2 normal approximation with the tie "
        "correction and the continuity correction, 3 the same without the continuity "
        "correction; alternative -1 less, 0 two-sided (twice the tail at max(U1, U2), "
        "at most 1), 1 greater"
    ),
    inputs=("x", "y", "method", "alternative"),
    outputs=("statistic", "pvalue"),
    tolerance=1e-13,
    tolerance_reason=(
        "For both outputs, by the rule of ttest_1samp. U1 is a sum of half-integer ranks, "
        "exact in double arithmetic, so its condition number in the data is 0 and only the "
        "few rounded operations in s move z, counted generously as n = n1 + n2; the exact "
        "p-value is a ratio of exact counts (kappa_p = 1), the asymptotic one adds "
        "kappa_p = |z phi(z) / Q| (mannwhitneyu_rounding_bound). The largest bound "
        "max(kappa_p, 1) n 2^-53 on this grid is 8.6e-14 (n1 = 10, n2 = 25, x above y, "
        "without continuity correction: z = 4.6); 1e-13 is the smallest round value above "
        "it (tested). SciPy's worst error is 2.1e-15"
    ),
    grid_description=(
        "x the first n1 and y the next n2 of v_i = (7919 i) mod 1000, x shifted by d: sizes "
        "(n1, n2) in {(3, 6), (6, 3), (4, 5), (5, 4), (8, 8), (8, 20), (20, 8)} with methods "
        "0 to 3, and {(9, 9), (10, 25), (25, 10)} with methods 0, 2 and 3, for d in "
        "{-1000.5, -200.5, 0.5, 150.5, 1000.5} (all values distinct); tied samples with "
        "v_i rounded down to a multiple of 200, sizes {(4, 5), (8, 8), (10, 25)}, d in "
        "{-1000, -200, 0, 200, 1000}, methods 0, 2 and 3; each with alternative -1 (less), "
        "0 (two-sided) and 1 (greater). U1 from 0 to 250 (54 points with a half-integer U1), p-value from "
        "3.2e-7 to 1"
    ),
    grid=_mwu_grid,
    mpmath_fn=_mp_mannwhitneyu,
    scipy_name="scipy.stats.mannwhitneyu",
    scipy_fn=_sc_mannwhitneyu,
)

# --- hypothesis tests: Fisher's exact test -------------------------------------

# 2x2 tables [[a, b], [c, d]] by their margins: row sums n1 = a + b and n2 = c + d,
# first column n = a + c, and a at the ends of its support, next to them, and at
# mode + k sd for k in _FISHER_SDS (sd of the hypergeometric distribution of a).
# Every row and column sum is positive (statrs and SciPy return p = 1 and odds
# ratio NaN otherwise, a convention, not a value), and b c > 0: an odds ratio of
# +inf is a division by zero, not a large finite value, yet the generator would tag
# it overflow, and at an overflow point a finite p-value, even a right one, is
# judged crate wrong.
# Above about 1030 values C(n1 + n2, n) exceeds the largest double, so the
# largest margins test that.
_FISHER_MARGINS = (
    (3, 5, 4),
    (5, 5, 5),
    (10, 10, 8),
    (8, 30, 12),
    (20, 20, 20),
    (50, 40, 30),
    (100, 100, 100),
    (200, 50, 60),
    (400, 400, 400),
    (600, 600, 300),
    (1000, 900, 700),
)
_FISHER_SDS = (-8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8)
# Published tables: SciPy's documentation (8, 2, 1, 5) and (6, 2, 1, 4), Fisher's
# tea tasting (3, 1, 1, 3), and the 7 tables of statrs's fisher test cases.
_FISHER_PUBLISHED = (
    (8, 2, 1, 5),
    (6, 2, 1, 4),
    (3, 1, 1, 3),
    (3, 5, 4, 50),
    (61, 118, 2, 1),
    (172, 46, 90, 127),
    (127, 38, 112, 43),
    (186, 177, 111, 154),
    (137, 49, 135, 183),
    (37, 115, 37, 152),
)


def _fisher_support(n1: int, n2: int, n: int) -> tuple[int, int]:
    return max(0, n - n2), min(n, n1)


def _fisher_near_ties(count: int = 12) -> list[tuple[int, int, int, int]]:
    """The first `count` tables (by total, then n1 <= n <= n2) whose probability
    is below another table's by at most 1e-4 relative: the two-sided p-value must
    leave that more probable table out (statrs's tie window is 1e-4 wide)."""
    found = []
    for total in itertools.count(4):
        for n1 in range(1, total // 2 + 1):
            for n in range(n1, total - n1 + 1):
                n2 = total - n1
                lo, hi = _fisher_support(n1, n2, n)
                w = [math.comb(n1, x) * math.comb(n2, n - x) for x in range(lo, hi + 1)]
                s = sorted(set(w))
                for u, v in itertools.pairwise(s):
                    a = lo + w.index(u)
                    table = (a, n1 - a, n - a, n2 - n + a)
                    if (v - u) * 10**4 <= v and table[1] * table[2] > 0:
                        found.append(table)
                        if len(found) == count:
                            return found
    raise AssertionError("unreachable")


def _fisher_tables() -> list[tuple[int, int, int, int]]:
    tables = list(_FISHER_PUBLISHED) + _fisher_near_ties()
    for n1, n2, n in _FISHER_MARGINS:
        lo, hi = _fisher_support(n1, n2, n)
        N = n1 + n2
        mean = n * n1 / N
        sd = math.sqrt(n * n1 * n2 * (N - n) / (N * N * (N - 1)))
        mode = (n + 1) * (n1 + 1) // (N + 2)
        picks = {lo, lo + 1, hi - 1, hi, mode}
        picks |= {round(mean + k * sd) for k in _FISHER_SDS}
        for a in sorted(v for v in picks if lo <= v <= hi):
            table = (a, n1 - a, n - a, n2 - (n - a))
            if table[1] * table[2] > 0 and table not in tables:
                tables.append(table)
    return tables


def _fisher_grid() -> list[tuple[float, float, float, float, float]]:
    return [(*map(float, t), alt) for t in _fisher_tables() for alt in _TTEST_ALTERNATIVES]


def _fisher_weights(a: int, b: int, c: int, d: int) -> tuple[list[int], int, int]:
    """w[x - lo] = C(n1, x) C(n2, n - x) for x in the support of the first cell: the
    number of tables with these margins and first cell x (its hypergeometric
    probability times C(n1 + n2, n)); also lo and the total C(n1 + n2, n)."""
    n1, n2, n = a + b, c + d, a + c
    lo, hi = _fisher_support(n1, n2, n)
    w = [math.comb(n1, x) * math.comb(n2, n - x) for x in range(lo, hi + 1)]
    return w, lo, math.comb(n1 + n2, n)


def _fisher_pvalue_counts(a: int, b: int, c: int, d: int, alternative: float) -> tuple[int, int]:
    """The p-value as (count, total), exactly: P(X <= a) (less), P(X >= a) (greater),
    or the sum of P(x) over every x with P(x) <= P(a) (two-sided), for X the first
    cell under the null hypothesis (hypergeometric, margins fixed)."""
    w, lo, total = _fisher_weights(a, b, c, d)
    i = a - lo
    if alternative < 0:
        return sum(w[: i + 1]), total
    if alternative > 0:
        return sum(w[i:]), total
    return sum(v for v in w if v <= w[i]), total


def _mp_fisher_exact(a, b, c, d, alternative):
    a, b, c, d = (int(v) for v in (a, b, c, d))
    odds = mpmath.mpf(a * d) / (b * c)  # the sample odds ratio; b c > 0 on the grid
    count, total = _fisher_pvalue_counts(a, b, c, d, float(alternative))
    return odds, min(mpmath.mpf(count) / total, mpmath.mpf(1))


def fisher_exact_rounding_bound(
    a: float, b: float, c: float, d: float, alternative: float
) -> float:
    """max(kappa_p, 1) (kappa + n) 2^-53 with kappa = 0 (whole-number inputs) and
    kappa_p = 1 (a ratio of exact counts); n = the number of hypergeometric terms
    in the support, a generous count of the rounded operations of a careful sum
    of the smaller tail. The odds ratio is a d / (b c), within 2 ulp."""
    lo, hi = _fisher_support(int(a + b), int(c + d), int(a + c))
    return float((hi - lo + 1) * mpmath.mpf(2) ** -53)


def _sc_fisher_exact(a: float, b: float, c: float, d: float, alternative: float):
    r = st.fisher_exact(
        [[int(a), int(b)], [int(c), int(d)]], alternative=_TTEST_ALTERNATIVES[alternative]
    )
    return float(r.statistic), float(r.pvalue)


FISHER_EXACT = Spec(
    name="fisher_exact",
    description=(
        "Fisher's exact test on the 2x2 table [[a, b], [c, d]]: the sample odds ratio "
        "a d / (b c) and the p-value from the hypergeometric "
        "distribution of the first cell with all margins fixed; alternative -1 less "
        "(P(X <= a)), 0 two-sided (the sum of P(x) over every x with P(x) <= P(a), "
        "compared exactly, at most 1), 1 greater (P(X >= a))"
    ),
    inputs=("a", "b", "c", "d", "alternative"),
    outputs=("statistic", "pvalue"),
    tolerance=1e-13,
    tolerance_reason=(
        "For both outputs, by the rule of ttest_1samp. The inputs are whole numbers "
        "(condition number 0) and the p-value is a ratio of exact counts (kappa_p = 1), so "
        "only rounding counts: n 2^-53 with n the number of hypergeometric terms in the "
        "support, a generous count for a careful sum of the smaller tail "
        "(fisher_exact_rounding_bound). The largest bound on this grid is 7.8e-14 "
        "(margins 1000, 900, 700: 701 terms); 1e-13 is the smallest round value above it "
        "(tested). The odds ratio is two products and a quotient, exact or within 2 ulp"
    ),
    grid_description=(
        "Published tables (SciPy's documentation, Fisher's tea tasting, statrs's "
        "fisher test cases: 10 tables), the first 12 tables (by total, then n1 <= n <= n2) whose "
        "probability is below another table's by at most 1e-4 relative (the first is "
        "(6, 9, 16, 41)), and tables by their margins (n1 = a + b, n2 = c + d, "
        "n = a + c) in {(3, 5, 4), (5, 5, 5), (10, 10, 8), (8, 30, 12), (20, 20, 20), "
        "(50, 40, 30), (100, 100, 100), (200, 50, 60), (400, 400, 400), (600, 600, 300), "
        "(1000, 900, 700)} with a at both ends of its support, next to them, at the mode, "
        "and at the mean + k sd for k in {+-8, +-4, +-2, +-1, +-0.5, 0}; each with "
        "alternative -1 (less), 0 (two-sided) and 1 (greater). No row or column sum is "
        "zero, and b c > 0 (an infinite odds ratio would be tagged overflow, and a "
        "finite p-value at an overflow point is judged crate wrong), which leaves out a at the top of its support"
    ),
    grid=_fisher_grid,
    mpmath_fn=_mp_fisher_exact,
    scipy_name="scipy.stats.fisher_exact",
    scipy_fn=_sc_fisher_exact,
)

# --- hypothesis tests: one-sample Kolmogorov-Smirnov test ----------------------

# Samples in (0, 1) tested against the uniform distribution on (0, 1), whose cdf
# F(x) = x is exact in statrs and SciPy, so only the statistic and the p-value are
# judged. z_i = ((7919 (i + 1)) mod 2048 + 1/2) / 2048, i = 0..n-1, are distinct
# exact doubles (7919 is odd); x_i = c z_i pushes the sample towards 0 (large D+)
# and x_i = 1 - c (1 - z_i) towards 1 (large D-), both exact for these c; the
# evenly spread x_i = (2 i - 1) / (2 n) give the smallest D.
# method: 0 less (statistic D- = max(F(x_(i)) - (i - 1) / n)), 1 greater
# (D+ = max(i / n - F(x_(i)))), both with the exact one-sided p-value; two-sided with
# D = max(D+, D-) and the p-value 2 exact (Marsaglia, Tsang and Wang 2003),
# 3 asymptotic (Kolmogorov's limit at D sqrt(n)), 4 approximate (twice the
# one-sided p-value, at most 1). statrs refuses the exact method from n = 170; the
# large sample (n = 1100) has the other four methods.
_KS_N = (1, 2, 3, 5, 10, 20, 50, 100, 150)
_KS_LARGE = ((1100, 1.0),)
_KS_SCALES = (1.0, 31 / 32, 7 / 8, 3 / 4, 5 / 8, 1 / 2)
_KS_METHODS = {
    0.0: ("less", "exact"),
    1.0: ("greater", "exact"),
    2.0: ("two-sided", "exact"),
    3.0: ("two-sided", "asymp"),
    4.0: ("two-sided", "approx"),
}


def _ks_sample(n: int, c: float, high: bool) -> tuple[float, ...]:
    z = [((7919 * (i + 1)) % 2048 + 0.5) / 2048 for i in range(n)]
    return tuple(1 - c * (1 - v) if high else c * v for v in z)


def _ks_spread_sample(n: int) -> tuple[float, ...]:
    """x_i = (2 i - 1) / (2 n), i = 1..n: as even as a sample can be, D = 1 / (2 n)
    (the smallest possible, up to the rounding of x_i), so D sqrt(n) is small."""
    return tuple((2 * i - 1) / (2 * n) for i in range(1, n + 1))


def _ks_grid() -> list[tuple[tuple[float, ...], float]]:
    grid = []
    shapes = [(1.0, False)] + [(c, high) for c in _KS_SCALES[1:] for high in (False, True)]
    for n in _KS_N:
        for c, high in shapes:
            x = _ks_sample(n, c, high)
            grid.extend((x, m) for m in _KS_METHODS)
        grid.extend((_ks_spread_sample(n), m) for m in _KS_METHODS)
    for n, c in _KS_LARGE:
        x = _ks_sample(n, c, False)
        grid.extend((x, m) for m in (0.0, 1.0, 3.0, 4.0))
    return grid


def _ks_d(x: list) -> tuple[Fraction, Fraction]:
    """(D+, D-) exactly, from the sorted sample as doubles: D+ = max(i / n - x_(i)),
    D- = max(x_(i) - (i - 1) / n), i = 1..n (F(x) = x)."""
    xs = sorted(Fraction(float(v)) for v in x)
    n = len(xs)
    d_plus = max(Fraction(i + 1, n) - v for i, v in enumerate(xs))
    d_minus = max(v - Fraction(i, n) for i, v in enumerate(xs))
    return d_plus, d_minus


def _ks_statistic(x: list, method: float) -> Fraction:
    d_plus, d_minus = _ks_d(x)
    if method == 0:
        return d_minus
    if method == 1:
        return d_plus
    return max(d_plus, d_minus)


def _mp_frac(q: Fraction) -> mpmath.mpf:
    return mpmath.mpf(q.numerator) / q.denominator


def _smirnov(n: int, d: Fraction) -> Fraction:
    """P(D+ >= d) for a sample of n, exactly (Birnbaum and Tingey 1951, eq. 3):
    d sum_{j=0}^{floor(n (1 - d))} C(n, j) (d + j / n)^(j - 1) (1 - d - j / n)^(n - j);
    0 for d > 1, 1 for d <= 0. With L the common denominator of d and 1 / n, every
    a_j = (d + j / n) L and b_j = (1 - d - j / n) L is a whole number, and the sum is
    (b_0^n + d L sum_{j>=1} C(n, j) a_j^(j - 1) b_j^(n - j)) / L^n (the j = 0 term
    times d is (1 - d)^n). Exact, so a p-value halfway between two doubles rounds
    the same way at every precision."""
    if d <= 0:
        return Fraction(1)
    if d > 1:
        return Fraction(0)
    big_l = d.denominator * n // math.gcd(d.denominator, n)
    dl = d.numerator * (big_l // d.denominator)
    step = big_l // n
    total = (big_l - dl) ** n
    for j in range(1, math.floor(n * (1 - d)) + 1):
        a = dl + j * step
        b = big_l - dl - j * step
        total += dl * math.comb(n, j) * a ** (j - 1) * b ** (n - j)
    return Fraction(total, big_l**n)


def _mp_smirnov(n: int, d: Fraction) -> mpmath.mpf:
    """P(D+ >= d) at the working precision (see _smirnov)."""
    return _mp_frac(_smirnov(n, d))


def _mp_kolmogorov_cdf_mtw(n: int, d: Fraction) -> mpmath.mpf:
    """P(D < d) for a sample of n by Marsaglia, Tsang and Wang (2003, section 2):
    n! / n^n (H^n)_{kk} with k = floor(n d) + 1, h = k - n d and the m x m matrix H,
    m = 2 k - 1. H has no negative entries, so H^n loses nothing; the vector
    H^n e_k is formed by n products."""
    k = math.floor(n * d) + 1
    m = 2 * k - 1
    h = Fraction(k) - n * d
    hm = _mp_frac(h)
    fact = [mpmath.factorial(i) for i in range(m + 1)]
    H = [
        [1 / fact[i - j + 1] if i - j + 1 >= 0 else mpmath.mpf(0) for j in range(m)]
        for i in range(m)
    ]
    for i in range(m):
        H[i][0] -= hm ** (i + 1) / fact[i + 1]
        H[m - 1][i] -= hm ** (m - i) / fact[m - i]
    if 2 * h > 1:
        H[m - 1][0] += _mp_frac(2 * h - 1) ** m / fact[m]
    v = [mpmath.mpf(0)] * m
    v[k - 1] = mpmath.mpf(1)
    for _ in range(n):
        # H[i][j] = 0 for j > i + 1 (lower Hessenberg).
        v = [mpmath.fsum(H[i][j] * v[j] for j in range(min(i + 2, m))) for i in range(m)]
    return v[k - 1] * mpmath.factorial(n) / mpmath.mpf(n) ** n


def _mp_kolmogorov_sf(n: int, d: Fraction) -> mpmath.mpf:
    """P(D >= d) for D = max(D+, D-): 2 P(D+ >= d) for d >= 1/2, where the two
    one-sided events are disjoint (up to a null set); otherwise 1 - P(D < d) with
    the digits that 1 - P loses added first (-log10 of the tail is at most about
    2 n d^2 log10(e) + 2)."""
    if d >= Fraction(1, 2):
        return 2 * _mp_smirnov(n, d)
    if 2 * n * d <= 1:
        return mpmath.mpf(1)  # D >= 1 / (2 n) for every sample
    extra = 20 + math.ceil(0.87 * n * d * d)
    with mpmath.extradps(extra):
        p = 1 - _mp_kolmogorov_cdf_mtw(n, d)
    return +p


def _mp_kolmogorov_limit_sf(t: mpmath.mpf) -> mpmath.mpf:
    """1 - K(t) for Kolmogorov's limit distribution K: 2 sum_{k>=1} (-1)^(k-1)
    e^(-2 k^2 t^2) for t >= 1, else 1 - sqrt(2 pi) / t sum_{k>=1}
    e^(-(2 k - 1)^2 pi^2 / (8 t^2)) (the Jacobi theta form; Feller 1948, eq. 2.6);
    both series converge fast in their range."""
    eps = mpmath.eps
    total = mpmath.mpf(0)
    k = 1
    if t >= 1:
        while True:
            term = mpmath.exp(-2 * k * k * t * t)
            total += term if k % 2 else -term
            if term < eps * total:
                return 2 * total
            k += 1
    while True:
        term = mpmath.exp(-((2 * k - 1) ** 2) * mpmath.pi**2 / (8 * t * t))
        total += term
        if term < eps * total:
            return 1 - mpmath.sqrt(2 * mpmath.pi) / t * total
        k += 1


def _mp_ks_pvalue(n: int, d: Fraction, method: float) -> mpmath.mpf:
    if method in (0, 1):
        return _mp_smirnov(n, d)
    if method == 2:
        return _mp_kolmogorov_sf(n, d)
    if method == 3:
        return _mp_kolmogorov_limit_sf(_mp_frac(d) * mpmath.sqrt(n))
    return min(2 * _mp_smirnov(n, d), mpmath.mpf(1))


def _mp_ks_1samp(x: list, method: mpmath.mpf):
    xs = [float(v) for v in x]  # the inputs are doubles
    d = _ks_statistic(xs, float(method))
    return _mp_frac(d), _mp_ks_pvalue(len(xs), d, float(method))


def ks_1samp_rounding_bound(x: tuple, method: float) -> float:
    """max(kappa_p, 1) (kappa_D + n) 2^-53. D = i / n - x_(i) (or x_(i) - (i - 1) / n)
    has one rounded quotient and one rounded difference, so it is off by up to
    (i / n + D) 2^-53 at its maximising i: kappa_D = (i / n + D) / D. kappa_p =
    |D p'(D) / p|, from a central difference at 40 digits; n counts generously the
    rounded operations of the p-value."""
    xs = sorted(float(v) for v in x)
    n = len(xs)
    d = _ks_statistic(xs, method)
    d_plus = _ks_d(xs)[0]
    if d == d_plus and method != 0:
        frac = max(
            Fraction(i + 1, n) for i, v in enumerate(xs) if Fraction(i + 1, n) - Fraction(v) == d
        )
    else:
        frac = max(Fraction(i, n) for i, v in enumerate(xs) if Fraction(v) - Fraction(i, n) == d)
    kappa_d = _mp_frac((frac + d) / d)
    with mpmath.workdps(40):
        step = d / 10**12
        p = _mp_ks_pvalue(n, d, method)
        p_up = _mp_ks_pvalue(n, d + step, method)
        p_down = _mp_ks_pvalue(n, d - step, method)
        kappa_p = abs((p_up - p_down) / (2 * _mp_frac(step)) * _mp_frac(d) / p) if p else 0
    return float(max(kappa_p, 1) * (kappa_d + n) * mpmath.mpf(2) ** -53)


def _sc_ks_1samp(x: list[float], method: float) -> tuple[float, float]:
    alternative, mode = _KS_METHODS[method]
    r = st.ks_1samp(x, st.uniform.cdf, alternative=alternative, method=mode)
    return float(r.statistic), float(r.pvalue)


KS_1SAMP = Spec(
    name="ks_1samp",
    description=(
        "One-sample Kolmogorov-Smirnov test of x against the uniform distribution on (0, 1) "
        "(F(x) = x): the statistic and its p-value; method 0 less (D- = max(F(x_(i)) - "
        "(i - 1) / n), p = P(D- >= d)), 1 greater (D+ = max(i / n - F(x_(i))), p = P(D+ >= d)), "
        "both exact (Birnbaum and Tingey 1951); two-sided with D = max(D+, D-) and p 2 exact "
        "P(D >= d) (Marsaglia, Tsang and Wang 2003), 3 asymptotic (1 - K(d sqrt(n)), K "
        "Kolmogorov's limit distribution), 4 approximate (min(2 P(D+ >= d), 1))"
    ),
    inputs=("x", "method"),
    outputs=("statistic", "pvalue"),
    tolerance=5e-12,
    tolerance_reason=(
        "For both outputs, by the rule of ttest_1samp. F is exact, and D = i / n - x_(i) (or "
        "x_(i) - (i - 1) / n) has one rounded quotient and one rounded difference, so "
        "kappa_D = (i / n + D) / D at the maximising i; the p-value adds kappa_p = "
        "|D p'(D) / p|, and n counts generously its rounded operations "
        "(ks_1samp_rounding_bound). The largest bound max(kappa_p, 1) (kappa_D + n) 2^-53 "
        "on this grid is 3.3e-12 (n = 150, x pushed towards 0 by c = 1/2, D = 0.53, "
        "p = 2.8e-40 one-sided, kappa_p = 197); 5e-12 is the smallest round value above it "
        "(tested). SciPy is wrong at 4 two-sided exact points at n = 150 (up to 2.6e-6: it "
        "switches to an approximation above n = 140); its worst below that is 3.6e-12"
    ),
    grid_description=(
        "z_i = ((7919 (i + 1)) mod 2048 + 1/2) / 2048, i = 0..n-1 (distinct exact doubles), "
        "and x = z, x = c z or x = 1 - c (1 - z) for c in {31/32, 7/8, 3/4, 5/8, 1/2} "
        "and the evenly spread x_i = (2 i - 1) / (2 n) (12 samples per n), for n in {1, 2, 3, "
        "5, 10, 20, 50, 100, 150} with methods 0 to 4, and x = z for n = 1100 with methods 0, "
        "1, 3 and 4 (statrs's exact method stops at n = 169). D from 3.7e-4 to 0.93, p-value "
        "from 2.8e-40 to 1"
    ),
    grid=_ks_grid,
    mpmath_fn=_mp_ks_1samp,
    scipy_name="scipy.stats.ks_1samp",
    scipy_fn=_sc_ks_1samp,
)

# --- betaincinv ----------------------------------------------------------------

# The betainc parameters; probabilities as for gamma_ppf, with fewer in the bulk.
_BETAINCINV_P = (0.0, 5e-324, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.25, 0.5,
                 0.75, 0.9, 0.99, 0.999, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1.0)  # fmt: skip


def _betaincinv_grid() -> list[tuple[float, float, float]]:
    pairs = itertools.product(_BETAINC_AB, repeat=2)
    return [(a, b, p) for a, b in pairs for p in _BETAINCINV_P]


def _mp_betainc_tail(a: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf, y: mpmath.mpf, lower: bool):
    """I_x(a, b) (lower) or 1 - I_x(a, b) (upper), with y = 1 - x given separately.

    The continued fraction is evaluated on whichever side of the switch point
    (a+1)/(a+b+2) it converges, at x or at y = 1 - x, so neither tail forms
    1 minus a number near 1 and both keep full relative precision.
    """
    if x <= (a + 1) / (a + b + 2):
        direct = _mp_betainc_cf(a, b, x)
        return direct if lower else 1 - direct
    mirrored = _mp_betainc_cf(b, a, y)  # = 1 - I_x(a, b)
    return 1 - mirrored if lower else mirrored


def _mp_betaincinv(a: mpmath.mpf, b: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with I_x(a, b) = p (0 <= p <= 1): x = 1 / (1 + e^-s) with s from
    _mp_betaincinv_logit; x = 0 at p = 0 and x = 1 at p = 1 (the ends of the support)."""
    if p == 0 or p == 1:
        return p
    s = _mp_betaincinv_logit(a, b, p)
    with mpmath.extradps(20):
        return 1 / (1 + mpmath.exp(-s))


def _mp_betaincinv_logit(a: mpmath.mpf, b: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """s = log(x / (1 - x)) for the x with I_x(a, b) = p (0 < p < 1), by
    Newton's method on a logarithm.

    In s = log(x / (1 - x)) both x = 1 / (1 + e^-s) and 1 - x = 1 / (1 + e^s)
    keep full relative precision, so quantiles near 0 and near 1 are equally
    accurate. For p <= 1/2 the root of log I = log p, for p > 1/2 the root of
    log(1 - I) = log(1 - p) (1 - p is exact for a double p at any working
    precision above 53 bits); both are smooth and monotone in s, with slope
    x^a (1-x)^b / (B(a, b) I) (respectively over 1 - I). Bracketing from the
    mean, bisection to within 1e-3 and then Newton's method, as for t_ppf
    and gamma_ppf.
    """
    lower = p <= mpmath.mpf(1) / 2
    q = p if lower else 1 - p
    sign = 1 if lower else -1
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    with mpmath.extradps(20):
        log_q = mpmath.log(q)
        log_beta = mpmath.log(mpmath.beta(a, b))

        def tail(s: mpmath.mpf) -> mpmath.mpf:
            x, y = 1 / (1 + mpmath.exp(-s)), 1 / (1 + mpmath.exp(s))
            return x, y, _mp_betainc_tail(a, b, x, y, lower)

        def g(s: mpmath.mpf) -> mpmath.mpf:
            # Increasing in s for both branches after the sign.
            return sign * (mpmath.log(tail(s)[2]) - log_q)

        s0 = mpmath.log(a / b)  # the mean a / (a + b)
        lo = hi = s0
        step = mpmath.mpf(1) / 1000
        while g(lo) > 0:
            lo -= step
            step *= 2
        step = mpmath.mpf(1) / 1000
        while g(hi) < 0:
            hi += step
            step *= 2
        while hi - lo > mpmath.mpf(1) / 1000:
            mid = (lo + hi) / 2
            if g(mid) > 0:
                hi = mid
            else:
                lo = mid
        s = (lo + hi) / 2
        for _ in range(100):
            x, y, t = tail(s)
            slope = mpmath.exp(a * mpmath.log(x) + b * mpmath.log(y) - log_beta) / t
            step = sign * (mpmath.log(t) - log_q) / slope
            s -= step
            if abs(step) <= eps:
                return s
    raise ArithmeticError(f"betaincinv({a}, {b}, {p}) did not converge")


BETAINCINV = Spec(
    name="betaincinv",
    description=(
        "Inverse of the regularized incomplete beta function: the x with I_x(a, b) = p "
        "(the beta distribution's quantile function)"
    ),
    inputs=("a", "b", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "a, b and p are exact doubles and 1 - p is exact for p >= 1/2, but in the lower tail the "
        "quantile is x = exp(L) with L = log(p a B(a, b)) / a approximately; rounding L (or 1/a, "
        "which multiplies it) to a double changes x by |L| times 2^-53. The largest |L| at a "
        "nonzero, normal quantile below 1 on this grid is 700.0 (a = 1, b = 1e4, p = 1e-300, "
        "x = 1e-304), so up to 7.8e-14 is unavoidable; near x = 1 the same holds for 1 - x, which "
        "moves x by less. 2e-13 allows 2.6x that and matches t_ppf and gamma_ppf. SciPy's worst "
        "error at the points it gets within tolerance is 2.0e-13 (a = 1000, b = 1e4, p = 0.75); "
        "it fails 105 points: NaN for 29 points with p <= 1e-300, about 2.2e-308 (just below the "
        "smallest normal double) where the quantile is smaller (41 points), errors from 1.2e-12 "
        "to 911x at 33 points with p <= 1e-300 (e.g. 0.375 instead of 0.477 at a = 1000, "
        "b = 0.5, p = 5e-324), and 2.8e-13 and 2.1e-13 at a = 1000, b = 1e4, p = 0.999 and "
        "1 - 1e-5"
    ),
    grid_description=(
        "a, b in {0.1, 0.5, 1, 2, 5, 10, 100, 1000, 10000} x p in {0, 5e-324, 1e-300, 1e-100, "
        "1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, and the doubles nearest 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15, 1}; the quantile is 0 at p = 0 and 1 at p = 1 (the ends of the "
        "support, so not poles); for small a the lower-tail quantiles are below the smallest "
        "double, and for small b the upper-tail quantiles round to 1"
    ),
    grid=_betaincinv_grid,
    mpmath_fn=_mp_betaincinv,
    scipy_name="scipy.special.betaincinv",
    scipy_fn=sc.betaincinv,
)

# --- beta_cdf and beta_sf ----------------------------------------------------

# The betainc parameters; x across both tails, the bulk and the ends of the
# support, plus the mean + k standard deviations (narrow for large a and b).
_BETA_X = (0.0, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9,
           0.99, 0.999, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1.0)  # fmt: skip
_BETA_K = (-10, -3, -1, 0, 1, 3, 10)


def _beta_x_grid() -> list[tuple[float, float, float]]:
    points = set()
    for a, b in itertools.product(_BETAINC_AB, repeat=2):
        points.update((a, b, x) for x in _BETA_X)
        mean = a / (a + b)
        sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
        points.update((a, b, x) for k in _BETA_K if 0 < (x := mean + k * sd) < 1)
    return sorted(points)


def _mp_beta_cdf(a: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """I_x(a, b), the beta distribution's cdf at x (0 <= x <= 1).

    The continued fraction on whichever side of the switch point converges
    (_mp_betainc_tail), with y = 1 - x exact, so neither tail loses relative
    precision; 0 and 1 at the ends of the support. Tested against
    mpmath.betainc, the positive-term series and closed forms.
    """
    if x == 0 or x == 1:
        return x
    return _mp_betainc_tail(a, b, x, mpmath.fsub(1, x, exact=True), lower=True)


def _mp_beta_sf(a: mpmath.mpf, b: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - I_x(a, b), the survival function, with full relative precision in both tails."""
    if x == 0 or x == 1:
        return 1 - x
    return _mp_betainc_tail(a, b, x, mpmath.fsub(1, x, exact=True), lower=False)


def beta_prefix_log_magnitude(a: float, b: float, x: float) -> float:
    """The largest |log| among the factors of the prefix x^a (1-x)^b / B(a, b)
    when formed separately: a ln x, b ln(1 - x), ln Gamma(a), ln Gamma(b) and
    ln Gamma(a + b). Rounding each to a double loses this times 2^-53."""
    a, b, x = mpmath.mpf(a), mpmath.mpf(b), mpmath.mpf(x)
    y = mpmath.fsub(1, x, exact=True)
    terms = (a * mpmath.log(x), b * mpmath.log(y), mpmath.loggamma(a), mpmath.loggamma(b),
             mpmath.loggamma(a + b))  # fmt: skip
    return float(max(abs(t) for t in terms))


_BETA_GRID_DESCRIPTION = (
    "a, b in {0.1, 0.5, 1, 2, 5, 10, 100, 1000, 10000} x x in {0, 1e-300, 1e-100, 1e-20, "
    "1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, and the doubles nearest "
    "1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1}, plus x = a/(a+b) + k sd (rounded to a double) for "
    "k in {-10, -3, -1, 0, 1, 3, 10} where 0 < x < 1, sd the standard deviation of Beta(a, b)"
)

BETA_CDF = Spec(
    name="beta_cdf",
    description=(
        "Beta distribution cumulative distribution function I_x(a, b) on [0, 1], "
        "including both tails and the ends of the support"
    ),
    inputs=("a", "b", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "Same function as betainc, so betainc's tolerance (accepted by the owner, to be "
        "revisited). The prefix x^a (1-x)^b / B(a, b) of the smaller of I and 1 - I is a "
        "product of factors whose logarithms (a ln x, b ln(1 - x), ln Gamma(a), ln Gamma(b), "
        "ln Gamma(a + b)) reach 1.78e5 in magnitude at a normal result on this grid "
        "(ln Gamma(2e4) at a = b = 1e4); rounding each to a double separately loses up to "
        "1.78e5 * 2^-53 = 2.0e-11, while combining them before rounding (Boost's Lanczos-based "
        "prefix) does not, so 1e-12 demands the careful method. SciPy's worst error here is "
        "1.1e-13 (a = b = 1000, x = 0.3), except 1.1e-9 at a = b = 1/2, x = 1 - 1e-15, where "
        "the closed form (2/pi) asin(sqrt x) confirms mpmath. The ends of the support give 0 "
        "and 1 exactly"
    ),
    grid_description=_BETA_GRID_DESCRIPTION,
    grid=_beta_x_grid,
    mpmath_fn=_mp_beta_cdf,
    scipy_name="scipy.stats.beta.cdf",
    scipy_fn=lambda a, b, x: st.beta.cdf(x, a, b),
)

BETA_SF = Spec(
    name="beta_sf",
    description="Beta distribution survival function 1 - I_x(a, b) on [0, 1]",
    inputs=("a", "b", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "1 - I, with beta_cdf's tolerance for the same reasons (largest prefix logarithm "
        "1.78e5 at a = b = 1e4). Computed as the complement 1 - I only where I is not close to "
        "1, the upper tail keeps full relative precision, which this tolerance demands down to "
        "the smallest normal double. SciPy's worst error here is 1.0e-15, except at a = b = 1/2, "
        "x = 1e-20, where it returns 1 for 1 - (2/pi) 1e-10 (6.4e-11 off; the closed form "
        "confirms mpmath)"
    ),
    grid_description="the beta_cdf grid",
    grid=_beta_x_grid,
    mpmath_fn=_mp_beta_sf,
    scipy_name="scipy.stats.beta.sf",
    scipy_fn=lambda a, b, x: st.beta.sf(x, a, b),
)

# --- F distribution (Fisher-Snedecor) -----------------------------------------

# Degrees of freedom d1, d2 (every pair): heavy tails (0.5, 1, 2), small and
# moderate integers, and large values where F concentrates near 1. With
# a = d1/2, b = d2/2 these are shapes the beta grid does not have.
_F_DF = (0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 1000.0, 1e5)
_F_X = (0.0, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0,
        100.0, 1e3, 1e5, 1e10, 1e20, 1e100, 1e300)  # fmt: skip
_F_K = (-10, -3, -1, 1, 3, 10)


def _f_x_grid() -> list[tuple[float, float, float]]:
    points = set()
    for d1, d2 in itertools.product(_F_DF, repeat=2):
        points.update((d1, d2, x) for x in _F_X)
        # Around 1 in steps of about the standard deviation of log F.
        sd = math.sqrt(2 / d1 + 2 / d2)
        points.update((d1, d2, x) for k in _F_K if (x := 1 + k * sd) > 0)
    return sorted(points)


def _mp_f_beta_args(d1: mpmath.mpf, d2: mpmath.mpf, x: mpmath.mpf):
    """a = d1/2, b = d2/2, z = d1 x / (d1 x + d2) and y = 1 - z = d2 / (d1 x + d2):
    F(x; d1, d2) = I_z(a, b) (A&S 26.6.2, DLMF 8.17 with the substitution
    z = d1 x / (d1 x + d2)). z and y are formed separately, so neither loses
    relative precision when the other is close to 1."""
    t = d1 * x
    return d1 / 2, d2 / 2, t / (t + d2), d2 / (t + d2)


def _mp_f_cdf(d1: mpmath.mpf, d2: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """F(x; d1, d2) for x >= 0, with full relative precision in both tails."""
    if x == 0:
        return mpmath.mpf(0)
    return _mp_betainc_tail(*_mp_f_beta_args(d1, d2, x), lower=True)


def _mp_f_sf(d1: mpmath.mpf, d2: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """1 - F(x; d1, d2) for x >= 0, with full relative precision in both tails."""
    if x == 0:
        return mpmath.mpf(1)
    return _mp_betainc_tail(*_mp_f_beta_args(d1, d2, x), lower=False)


def f_condition(d1: float, d2: float, x: float, upper: bool) -> mpmath.mpf:
    """|x f(x) / F(x)| (upper: over 1 - F(x)): the relative change of F per
    relative change of x, z or y. x f(x) = z^a y^b / B(a, b), so it is also the
    beta prefix over the tail."""
    d1, d2, x = mpmath.mpf(d1), mpmath.mpf(d2), mpmath.mpf(x)
    a, b, z, y = _mp_f_beta_args(d1, d2, x)
    prefix = mpmath.exp(a * mpmath.log(z) + b * mpmath.log(y) - mpmath.log(mpmath.beta(a, b)))
    tail = _mp_f_sf(d1, d2, x) if upper else _mp_f_cdf(d1, d2, x)
    return prefix / tail


def _f_ppf_grid() -> list[tuple[float, float, float]]:
    pairs = itertools.product(_F_DF, repeat=2)
    return [(d1, d2, p) for d1, d2 in pairs for p in _BETAINCINV_P]


def _f_ppf_tags(d1: float, d2: float, p: float) -> list[str]:
    # The quantile is +inf at p = 1; at p = 0 it is 0, the end of the support.
    return [POLE] if p == 1.0 else []


def _mp_f_ppf(d1: mpmath.mpf, d2: mpmath.mpf, p: mpmath.mpf) -> mpmath.mpf:
    """The x with F(x; d1, d2) = p (0 <= p < 1).

    With z = I^-1_p(d1/2, d2/2), x = (d2 / d1) z / (1 - z) = (d2 / d1) e^s for
    s = log(z / (1 - z)) from _mp_betaincinv_logit, so neither z nor 1 - z is
    rounded on the way and quantiles in both tails keep full relative precision.
    """
    if p == 0:
        return mpmath.mpf(0)
    s = _mp_betaincinv_logit(d1 / 2, d2 / 2, p)
    with mpmath.extradps(20):
        return d2 / d1 * mpmath.exp(s)


_F_DF_DESCRIPTION = "d1, d2 in {0.5, 1, 2, 5, 10, 30, 100, 1000, 1e5} (every pair)"

F_CDF = Spec(
    name="f_cdf",
    description=(
        "F (Fisher-Snedecor) distribution cumulative distribution function F(x; d1, d2) "
        "= I_z(d1/2, d2/2) with z = d1 x / (d1 x + d2)"
    ),
    inputs=("d1", "d2", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "Same function as betainc (F(x) = I_z(d1/2, d2/2)), so betainc's tolerance (accepted by "
        "the owner, to be revisited). Two error sources: (1) forming "
        "z = d1 x / (d1 x + d2) or y = 1 - z = d2 / (d1 x + d2) directly costs a few roundings, "
        "which change F by the condition number |x f(x) / F(x)| (= z^a y^b / (B(a, b) F)) times "
        "2^-53 each; its largest value at a normal result on this grid is 1648 (d1 = d2 = 1e5, "
        "x = 0.937), so about 1.8e-13 per rounding. Forming y as 1 - z instead loses "
        "2^-53 z / y relative to y, unbounded as x grows, and is not needed. (2) The prefix "
        "z^a y^b / B(a, b) has logarithms up to 1.05e6 (ln Gamma(1e5) at d1 = d2 = 1e5); rounding "
        "them separately costs up to 1.2e-10, combining them before rounding (Boost's "
        "Lanczos-based prefix) does not, so 1e-12 demands the careful method. SciPy's worst error "
        "at the points it gets within tolerance is 9.9e-13 (d1 = 1e5, d2 = 10, x = 0.553); it "
        "fails 6: 1.1e-12 to 1.7e-12 at d1 = 1e5 with d2 = 10 or 30 (it forms z near 1 and "
        "betainc rounds 1 - z), and 1 instead of 1 - 6.4e-11 at d1 = d2 = 1, x = 1e20 (the "
        "closed form (2/pi) atan(sqrt x) confirms mpmath). x = 0 gives 0 exactly"
    ),
    grid_description=(
        f"{_F_DF_DESCRIPTION} x x in {{0, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, "
        "0.5, 1, 2, 5, 10, 100, 1e3, 1e5, 1e10, 1e20, 1e100, 1e300}, plus x = 1 + k s (rounded "
        "to a double) for k in {-10, -3, -1, 1, 3, 10} where x > 0, s = sqrt(2/d1 + 2/d2) "
        "(about the standard deviation of log F)"
    ),
    grid=_f_x_grid,
    mpmath_fn=_mp_f_cdf,
    scipy_name="scipy.stats.f.cdf",
    scipy_fn=lambda d1, d2, x: st.f.cdf(x, d1, d2),
)

F_SF = Spec(
    name="f_sf",
    description="F (Fisher-Snedecor) distribution survival function 1 - F(x; d1, d2)",
    inputs=("d1", "d2", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "1 - F, with f_cdf's tolerance for the same reasons: the condition number |x f / (1 - F)| "
        "reaches 1549 at a normal result (d1 = d2 = 1e5, x = 1.063), about 1.7e-13 per rounding "
        "of z or y, and the prefix logarithms reach 1.05e6, so separately rounded logarithms can "
        "cost 1.2e-10. Computed from y = d2 / (d1 x + d2) directly, the upper tail keeps full "
        "relative precision for every x, which this tolerance demands down to the smallest "
        "normal double. SciPy's worst error at the points it gets within tolerance is 8.6e-13 "
        "(d1 = 30, d2 = 1e5, x = 5); it fails 11: 1.0e-12 to 2.4e-12 at d2 = 1e5 with d1 = 10 or "
        "30, and 1 instead of 1 - 6.4e-11 at d1 = d2 = 1, x = 1e-20 (the closed form confirms "
        "mpmath)"
    ),
    grid_description="the f_cdf grid",
    grid=_f_x_grid,
    mpmath_fn=_mp_f_sf,
    scipy_name="scipy.stats.f.sf",
    scipy_fn=lambda d1, d2, x: st.f.sf(x, d1, d2),
)

F_PPF = Spec(
    name="f_ppf",
    description="F (Fisher-Snedecor) distribution quantile function F^-1(p; d1, d2)",
    inputs=("d1", "d2", "p"),
    tolerance=2e-13,
    tolerance_reason=(
        "d1, d2 and p are exact doubles and 1 - p is exact for p >= 1/2, but in the tails the "
        "quantile is x = exp(L) with L of the order of log(p) / a (lower tail) or "
        "-log(1 - p) / b (upper tail); rounding L (or 1/a, 1/b, d2/d1, which multiply it) to a "
        "double changes x by |L| times 2^-53. The largest |L| at a positive, normal quantile on "
        "this grid is 690.8 (d1 = 2, p = 1e-300, x = 1e-300), so up to 7.7e-14 is unavoidable; "
        "2e-13 allows 2.6x that and matches betaincinv. SciPy's worst error at the points it "
        "gets within tolerance is 1.8e-13 (d1 = d2 = 1000, p = 1e-100); it fails 101 points, "
        "all with p <= 1e-100: NaN at 27 points, numbers from 1.1e-308 to 4.5e-303 "
        "where the quantile is 0 or subnormal (41 points), and errors from 8.5e-12 to 911x at "
        "33 points (e.g. 3e-5 instead of 6.7e-4 at d1 = 1e5, d2 = 1, p = 5e-324)"
    ),
    grid_description=(
        f"{_F_DF_DESCRIPTION} x p in {{0, 5e-324, 1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 1e-3, "
        "0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, and the doubles nearest 1 - 1e-5, "
        "1 - 1e-10, 1 - 1e-15, 1}; the quantile is 0 at p = 0 (the end of the support) and "
        "p = 1 is tagged pole (+inf); quantiles beyond the largest double are tagged overflow"
    ),
    grid=_f_ppf_grid,
    mpmath_fn=_mp_f_ppf,
    scipy_name="scipy.stats.f.ppf",
    scipy_fn=lambda d1, d2, p: st.f.ppf(p, d1, d2),
    input_tags=_f_ppf_tags,
)

# --- gammaln and gammasgn --------------------------------------------------

# Near the zeros of log|Gamma| at 1 and 2: x - 1 and x - 2 are exact there.
_GAMMALN_NEAR_ZEROS = (1e-5, 1e-10, 1e-15)
# The overflow edge: log Gamma(x) exceeds the largest double above about 2.55e305.
_GAMMALN_HUGE = (1e305, 2e305, 2.5e305, 2.6e305, 1e306, 1e307, 1e308, 1.7976931348623157e308)
# Large negative non-integers (exact doubles), reached through the reflection formula.
_GAMMALN_NEGATIVE = (-100.5, -1000.25, -10000.5, -1000000.5, -999999999999999.5)


def _gammaln_grid() -> list[tuple[float]]:
    # The gamma grid (x = -30..180 step 0.25 and +-10^-k), plus decades up to
    # 1e300, the overflow edge, points next to the zeros at 1 and 2, and large
    # negative non-integers.
    points = {x for (x,) in _gamma_grid()}
    points.update(_decades())
    points.update(_GAMMALN_HUGE)
    points.update(_GAMMALN_NEGATIVE)
    for zero in (1.0, 2.0):
        for d in _GAMMALN_NEAR_ZEROS:
            points.update((zero + d, zero - d))
    return [(x,) for x in sorted(points)]


def _mp_gammaln(x: mpmath.mpf) -> mpmath.mpf:
    # log|Gamma(x)|; mpmath's loggamma is complex for negative x, with this
    # real part.
    return mpmath.re(mpmath.loggamma(x))


def gammaln_reflection_amplification(x: float) -> mpmath.mpf:
    """For x < 0: (|log(pi / |sin(pi x)|)| + |log Gamma(1 - x)|) / |log|Gamma(x)||.

    The reflection formula log|Gamma(x)| = log(pi / |sin(pi x)|) - log Gamma(1 - x)
    subtracts two terms; rounding each to a double costs this many times 2^-53
    relative to the result. It is large only near the zeros of log|Gamma|.
    """
    x = mpmath.mpf(x)
    t1 = mpmath.log(mpmath.pi / abs(mpmath.sinpi(x)))
    t2 = mpmath.loggamma(1 - x)
    return (abs(t1) + abs(t2)) / abs(_mp_gammaln(x))


GAMMALN = Spec(
    name="gammaln",
    description="Logarithm of the absolute value of the gamma function, log|gamma(x)|, for real x",
    inputs=("x",),
    tolerance=2e-13,
    tolerance_reason=(
        "For x > 0 the value is well conditioned once x - 1 and x - 2 are formed exactly "
        "(Sterbenz's lemma near the zeros at 1 and 2, where careful implementations such as "
        "fdlibm and Boost expand in x - 1 and x - 2; taking log(gamma(x)) there instead loses "
        "everything), and Stirling's series for large x costs a few ulps. For x < 0 the "
        "reflection formula log(pi / |sin(pi x)|) - log Gamma(1 - x) subtracts two terms; "
        "rounding them costs (|term 1| + |term 2|) / |result| times 2^-53, at most 664 on this "
        "grid (x = -2.75, near the zero at -2.7477), so up to 7.4e-14 is unavoidable. 2e-13 "
        "allows 2.7x that and matches the quantile tolerances. SciPy's worst error is 1.8e-15 "
        "except at 8 of the 12 points next to 1 and 2, where it loses digits as log(gamma(x)) "
        "would: 1.5e-12 at 1 +- 1e-5, up to 8.0e-7 at 1 +- 1e-10 and 2 - 1e-10, and 29% at "
        "the double nearest 2 - 1e-15 (the series DLMF 5.7.3 confirms mpmath)"
    ),
    grid_description=(
        "x = -30..180 step 0.25 (poles at 0, -1, ..., -30) and +-10^k for k = -300..-10 step "
        "10 (the gamma grid), plus 10^k for k = 0..300 step 10, x = 1e305, 2e305, 2.5e305, "
        "2.6e305, 1e306, 1e307, 1e308 and the largest double (log Gamma overflows above about "
        "2.55e305), 1 +- d and 2 +- d for d in {1e-5, 1e-10, 1e-15} (rounded to doubles), and "
        "x = -100.5, -1000.25, -10000.5, -1000000.5, -999999999999999.5"
    ),
    grid=_gammaln_grid,
    mpmath_fn=_mp_gammaln,
    scipy_name="scipy.special.gammaln",
    scipy_fn=sc.gammaln,
    input_tags=_gamma_tags,
)


def _gammasgn_grid() -> list[tuple[float]]:
    # The poles are left out: log|Gamma| is +inf there, and the conventions
    # accept an infinity of either sign at a pole, so no sign is wrong there.
    return [(x,) for (x,) in _gammaln_grid() if not _gamma_tags(x)]


def _mp_gammasgn(x: mpmath.mpf) -> mpmath.mpf:
    return mpmath.sign(mpmath.gamma(x))


GAMMASGN = Spec(
    name="gammasgn",
    description="Sign of the gamma function (+1 or -1) for real x, as returned with log|gamma(x)|",
    inputs=("x",),
    tolerance=0.0,
    tolerance_reason=(
        "The sign is exactly +1 or -1 at every input of the grid, so only the exact value is "
        "right. For x < 0 it is (-1)^ceil(-x), for x > 0 it is +1"
    ),
    grid_description=(
        "the gammaln grid without its poles (0, -1, ..., -30), where gamma has no sign and the "
        "conventions accept an infinity of either sign"
    ),
    grid=_gammasgn_grid,
    mpmath_fn=_mp_gammasgn,
    scipy_name="scipy.special.gammasgn",
    scipy_fn=sc.gammasgn,
)


# --- gammainc_large --------------------------------------------------------

# Shapes around and above 1e8; for 1e8 and 1e9 also the bulk x = a + k sqrt(a)
# (exact for a = 1e8; x - a is exact for every point by Sterbenz's lemma).
_GAMMAINC_LARGE_A = (1e8, 1e9, 1e10, 1e15, 1e20)
_GAMMAINC_LARGE_BULK_A = (1e8, 1e9)
_GAMMAINC_LARGE_K = (-30, -10, -3, -1, 0, 1, 3, 10, 30)


def _gammainc_large_grid() -> list[tuple[float, float]]:
    points = set()
    for a in _GAMMAINC_LARGE_A:
        points.update((a, x) for x in (1e8, 2e8, a / 2, 0.9 * a, 1.1 * a, 2 * a, 1e300))
    for a in _GAMMAINC_LARGE_BULK_A:
        points.update((a, a + k * math.sqrt(a)) for k in _GAMMAINC_LARGE_K)
    return sorted(points)


GAMMAINC_LARGE = Spec(
    name="gammainc_large",
    description=("Regularized lower incomplete gamma function P(a, x) for large shapes a >= 1e8"),
    inputs=("a", "x"),
    tolerance=1e-12,
    tolerance_reason=(
        "The function and method of gammainc and gamma_cdf, with their tolerance. The condition "
        "number |x f(x) / P| grows like k sqrt(a) at x = a + k sqrt(a) (9.5e5 at a = 1e9, "
        "k = -30), so a method that rounds x / a would lose about 1e-10; one that forms x - a "
        "exactly (Sterbenz's lemma) and expands in (x - a) / a (Temme, DLMF 8.12; Boost) loses "
        "about k^2 2^-53 (1e-13 at k = 30). Far from the bulk P rounds to 0 or 1. SciPy's worst "
        "error is 2.7e-15 except in the lower tail at x = a - 10 sqrt(a) and a - 30 sqrt(a): "
        "13% and 0.24% off at a = 1e8, 53% and 15% off at a = 1e9 (P about 7e-24 and 2e-198; "
        "numerical integration confirms mpmath)"
    ),
    grid_description=(
        "a in {1e8, 1e9, 1e10, 1e15, 1e20} x x in {1e8, 2e8, a/2, 0.9a, 1.1a, 2a, 1e300}, plus "
        "x = a + k sqrt(a) for a in {1e8, 1e9} and k in {-30, -10, -3, -1, 0, 1, 3, 10, 30}"
    ),
    grid=_gammainc_large_grid,
    mpmath_fn=_mp_gamma_cdf,
    scipy_name="scipy.special.gammainc",
    scipy_fn=sc.gammainc,
)

# --- digamma and trigamma ---------------------------------------------------

# The positive zero of digamma, x0 = 1.46163214496836234126...: the double
# nearest it, and x0 +- d for d in _DIGAMMA_NEAR_ROOT (rounded to doubles).
_DIGAMMA_ROOT = 1.4616321449683623
_DIGAMMA_NEAR_ROOT = (1e-5, 1e-10)


def _polygamma_grid() -> list[tuple[float]]:
    # The gamma grid (x = -30..180 step 0.25 and +-10^-k), plus decades up to
    # 1e300, the gammaln grid's huge values and large negative non-integers.
    points = {x for (x,) in _gamma_grid()}
    points.update(_decades())
    points.update(_GAMMALN_HUGE)
    points.update(_GAMMALN_NEGATIVE)
    return [(x,) for x in sorted(points)]


def _digamma_grid() -> list[tuple[float]]:
    points = {x for (x,) in _polygamma_grid()}
    points.add(_DIGAMMA_ROOT)
    for d in _DIGAMMA_NEAR_ROOT:
        points.update((_DIGAMMA_ROOT + d, _DIGAMMA_ROOT - d))
    return [(x,) for x in sorted(points)]


def digamma_reflection_amplification(x: float) -> mpmath.mpf:
    """For x < 0: (|psi(1 - x)| + |pi cot(pi x)|) / |psi(x)|.

    The reflection formula psi(x) = psi(1 - x) - pi cot(pi x) subtracts two
    terms; rounding each to a double costs this many times 2^-53 relative to
    the result. It is large only near the zeros of psi.
    """
    x = mpmath.mpf(x)
    t1 = mpmath.digamma(1 - x)
    t2 = mpmath.pi * mpmath.cospi(x) / mpmath.sinpi(x)
    return (abs(t1) + abs(t2)) / abs(mpmath.digamma(x))


DIGAMMA = Spec(
    name="digamma",
    description="Digamma function psi(x) = gamma'(x) / gamma(x) for real x",
    inputs=("x",),
    tolerance=4e-13,
    tolerance_reason=(
        "For x > 0 the value is well conditioned except next to the zero x0 = 1.4616..., where "
        "careful implementations (Boost, SciPy) expand in x - x0 with x0 held to extra "
        "precision, so the exact input still gives a result to a few ulps; the grid has the "
        "double nearest x0 and x0 +- 1e-5, 1e-10 to test this. For x < 0 the reflection formula "
        "psi(1 - x) - pi cot(pi x) subtracts two terms; rounding them costs (|term 1| + "
        "|term 2|) / |result| times 2^-53, at most 1313 on this grid (x = -22.75, near a zero "
        "at -22.7502), so up to 1.5e-13 is unavoidable. 4e-13 allows 2.7x that, as for "
        "gammaln. SciPy's worst error is 2.8e-13, at that point; elsewhere 2.2e-14"
    ),
    grid_description=(
        "x = -30..180 step 0.25 (poles at 0, -1, ..., -30) and +-10^k for k = -300..-10 step "
        "10 (the gamma grid), plus 10^k for k = 0..300 step 10, x = 1e305, 2e305, 2.5e305, "
        "2.6e305, 1e306, 1e307, 1e308 and the largest double, x = -100.5, -1000.25, -10000.5, "
        "-1000000.5, -999999999999999.5, and next to the positive zero x0 = 1.46163...: the "
        "double nearest x0 and x0 +- d for d in {1e-5, 1e-10} (rounded to doubles)"
    ),
    grid=_digamma_grid,
    mpmath_fn=mpmath.digamma,
    scipy_name="scipy.special.digamma",
    scipy_fn=sc.digamma,
    input_tags=_gamma_tags,
)


def _trigamma_grid() -> list[tuple[float]]:
    # SciPy's polygamma(1, x) sums one term per unit step up from negative x
    # (0.05 s at x = -1e6), so at x = -999999999999999.5 it does not return;
    # every point needs a SciPy value, so that one point is left out.
    return [(x,) for (x,) in _polygamma_grid() if x != -999999999999999.5]


def _mp_trigamma(x: mpmath.mpf) -> mpmath.mpf:
    # mpmath.psi(1, x) recurs one step at a time up from negative x (5 s at
    # x = -1e6, endless at -1e15), so negative x uses the reflection formula
    # psi'(x) = pi^2 / sin^2(pi x) - psi'(1 - x) (derivative of DLMF 5.5.4).
    if x < 0:
        return mpmath.pi**2 / mpmath.sinpi(x) ** 2 - mpmath.psi(1, 1 - x)
    return mpmath.psi(1, x)


def _sc_trigamma(x: float) -> float:
    return sc.polygamma(1, x)


TRIGAMMA = Spec(
    name="trigamma",
    description="Trigamma function psi'(x), the derivative of digamma, for real x",
    inputs=("x",),
    tolerance=1e-13,
    tolerance_reason=(
        "Trigamma is positive and has no zeros. For x > 0 the recurrence and the asymptotic "
        "series add terms of one sign. For x < 0 the reflection formula pi^2 / sin^2(pi x) - "
        "psi'(1 - x) subtracts at most psi'(1) = pi^2 / 6 from at least pi^2, so rounding is "
        "amplified at most 1.4 times, provided sin(pi x) is formed from the exact input (as "
        "sinpi does). 1e-13 is the tolerance of the other well-conditioned functions"
    ),
    grid_description=(
        "the digamma grid without the points next to its zero and without x = "
        "-999999999999999.5, where SciPy's polygamma(1, x) does not return (it sums one term "
        "per unit step from x up)"
    ),
    grid=_trigamma_grid,
    mpmath_fn=_mp_trigamma,
    scipy_name="scipy.special.polygamma(1, x)",
    scipy_fn=_sc_trigamma,
    input_tags=_gamma_tags,
)

ALL: tuple[Spec, ...] = (
    ERFCX, ERFC, GAMMA, GAMMAINC, BETAINC, BETAINC_LARGE, BETAINC_SWITCH, KV, YV,
    YV_LARGE, JV, JV_LARGE, IV, IV_LARGE, IV_NONINT, NORM_CDF, NORM_SF, NORM_PPF,
    T_CDF, T_SF, T_PPF, GAMMA_CDF, GAMMA_SF, GAMMA_PPF,
    INVGAMMA_CDF, INVGAMMA_SF, INVGAMMA_PPF, CHI_CDF, CHI_SF, CHI_PPF, EXPON_CDF, EXPON_SF, EXPON_PPF,
    WEIBULL_CDF, WEIBULL_SF, WEIBULL_PPF, LOGNORM_CDF, LOGNORM_SF, LOGNORM_PPF,
    PARETO_CDF, PARETO_SF, PARETO_PPF,
    GUMBEL_CDF, GUMBEL_SF, GUMBEL_PPF, LAPLACE_CDF, LAPLACE_SF, LAPLACE_PPF,
    CAUCHY_CDF, CAUCHY_SF, CAUCHY_PPF, LEVY_CDF, LEVY_SF, LEVY_PPF,
    BETAINCINV, BETA_CDF, BETA_SF, TTEST_1SAMP, CHISQUARE, F_ONEWAY,
    SKEWTEST, MANNWHITNEYU, FISHER_EXACT, KS_1SAMP,
    F_CDF, F_SF, F_PPF,
    GAMMALN, GAMMASGN, GAMMAINC_LARGE, DIGAMMA, TRIGAMMA,
)  # fmt: skip
