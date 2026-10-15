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
    INVGAMMA_CDF, INVGAMMA_SF, INVGAMMA_PPF, EXPON_CDF, EXPON_SF, EXPON_PPF,
    BETAINCINV, BETA_CDF, BETA_SF,
    F_CDF, F_SF, F_PPF,
    GAMMALN, GAMMASGN, GAMMAINC_LARGE, DIGAMMA, TRIGAMMA,
)  # fmt: skip
