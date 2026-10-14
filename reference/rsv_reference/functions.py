"""The registry of reference datasets. One entry per mathematical function."""

import itertools
import math

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
    """The x with P(a, x) = p (0 < p < 1), by Newton's method on a logarithm.

    For p <= 1/2 the root of log P(a, e^u) = log p, for p > 1/2 the root of
    log Q(a, e^u) = log(1 - p) (1 - p is exact for a double p at any working
    precision above 53 bits). In u = log x both sides are smooth and monotone,
    with slopes x f / P and -x f / Q. The search starts from the
    Wilson-Hilferty approximation (or, where that is not positive, from the
    lower-tail power law P = x^a / Gamma(a + 1)), brackets the root, bisects
    to within 1e-3 and then converges quadratically. ppf(0) = 0.
    """
    if p == 0:
        return mpmath.mpf(0)
    lower = p <= mpmath.mpf(1) / 2
    q = p if lower else 1 - p
    sign = 1 if lower else -1
    eps = mpmath.mpf(10) ** (-mpmath.mp.dps - 5)
    with mpmath.extradps(20):
        log_q = mpmath.log(q)

        def g(u: mpmath.mpf) -> mpmath.mpf:
            # Increasing in u for both branches after the sign.
            x = mpmath.exp(u)
            tail = _mp_gamma_cdf(a, x) if lower else _mp_gamma_sf(a, x)
            return sign * (mpmath.log(tail) - log_q)

        z = -mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * p)  # the normal quantile
        w = 1 - 1 / (9 * a) + z / (3 * mpmath.sqrt(a))
        u0 = 3 * mpmath.log(w) + mpmath.log(a) if w > 0 else (log_q + mpmath.loggamma(a + 1)) / a
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
    raise ArithmeticError(f"gamma_ppf({a}, {p}) did not converge")


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

ALL: tuple[Spec, ...] = (
    ERFCX, ERFC, GAMMA, GAMMAINC, BETAINC, BETAINC_LARGE, BETAINC_SWITCH, KV, YV,
    YV_LARGE, JV, JV_LARGE, IV, IV_LARGE, IV_NONINT, NORM_CDF, NORM_SF, NORM_PPF,
    T_CDF, T_SF, T_PPF, GAMMA_CDF, GAMMA_SF, GAMMA_PPF,
)  # fmt: skip
