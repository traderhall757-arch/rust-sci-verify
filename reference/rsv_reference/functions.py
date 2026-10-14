"""The registry of reference datasets. One entry per mathematical function."""

import itertools
import math

import mpmath
import scipy.special as sc

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

ALL: tuple[Spec, ...] = (
    ERFCX, ERFC, GAMMA, GAMMAINC, BETAINC, BETAINC_LARGE, BETAINC_SWITCH, KV, YV,
    YV_LARGE,
)  # fmt: skip
