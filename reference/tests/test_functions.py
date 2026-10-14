"""Each mpmath route used as the arbiter, checked against an independent one.

Agreement between 50 and 100 digits does not catch a wrong formula, so every
dataset's mpmath function is compared here with a different method (an exact
identity, a finite sum, or a different mpmath function) where both are valid.
"""

import math

import mpmath
import pytest

from rsv_reference import functions
from rsv_reference.spec import POLE


def close(a, b, digits=40):
    return abs(a - b) <= abs(b) * mpmath.mpf(10) ** -digits


def mp(*args):
    return [mpmath.mpf(a) for a in args]


def test_gamma_matches_factorials_and_half_integers():
    with mpmath.workdps(50):
        for n in (1, 2, 7, 30, 171):
            assert close(mpmath.gamma(n), mpmath.factorial(n - 1))
        for n in (0, 3, 25):
            # gamma(n + 1/2) = (2n)! sqrt(pi) / (4^n n!)
            exact = mpmath.factorial(2 * n) * mpmath.sqrt(mpmath.pi)
            exact /= mpmath.mpf(4) ** n * mpmath.factorial(n)
            assert close(mpmath.gamma(n + mpmath.mpf(0.5)), exact)


def test_gamma_poles_are_tagged_by_rule():
    assert functions._gamma_tags(-6.0) == [POLE]
    assert functions._gamma_tags(0.0) == [POLE]
    assert functions._gamma_tags(-6.5) == []
    assert functions._gamma_tags(6.0) == []


@pytest.mark.parametrize("x", [0.5, 3.0, 26.75, 1e3])
def test_erfc_matches_upper_incomplete_gamma(x):
    # erfc(x) = Gamma(1/2, x^2) / sqrt(pi) for x > 0 (DLMF 7.11.2).
    with mpmath.workdps(50):
        (x,) = mp(x)
        route = mpmath.gammainc(mpmath.mpf(0.5), x * x) / mpmath.sqrt(mpmath.pi)
        assert close(functions._mp_erfc(x), route)


def test_erfc_large_x_branch_matches_direct_formula():
    with mpmath.workdps(50):
        # Just above 1e6 the series branch is used; mpmath's erfc is still valid there.
        x = mpmath.mpf(10) ** 6 * (1 + mpmath.mpf(10) ** -45)
        # Both routes round x^2 (about 1e12) at 50 digits inside exp, which
        # costs about 1e-38; 1e-30 is still far below double precision.
        assert abs(functions._mp_erfc(x) / mpmath.erfc(x) - 1) < mpmath.mpf(10) ** -30
        assert functions._mp_erfc(mpmath.mpf(1e300)) > 0  # tiny, not an error


@pytest.mark.parametrize("a, x", [(1, 0.5), (5, 2.0), (500, 100.0), (1000, 2000.0), (20, 20.0)])
def test_gammainc_matches_finite_sum_for_integer_a(a, x):
    # For integer a, Q(a, x) = e^-x sum_{k<a} x^k / k! (DLMF 8.4.10). At 400
    # digits, 1 - Q keeps enough digits even where P is about 1e-178.
    with mpmath.workdps(400):
        xm = mpmath.mpf(x)
        q = mpmath.exp(-xm) * mpmath.fsum(xm**k / mpmath.factorial(k) for k in range(a))
        exact = 1 - q
    with mpmath.workdps(50):
        assert close(functions._mp_gammainc(*mp(a, x)), exact)


@pytest.mark.parametrize(
    "a, b, x",
    [(2, 3, 0.3), (1000, 1000, 0.5), (500, 500, 0.49), (10000, 10000, 0.5), (3, 10000, 0.01)],
)
def test_betainc_series_matches_binomial_sum_for_integer_a_b(a, b, x):
    # For integer a, b: I_x(a, b) = sum_{j=a}^{n} C(n, j) x^j (1-x)^(n-j), n = a+b-1
    # (DLMF 8.17.5). Positive terms, evaluated at 60 digits.
    with mpmath.workdps(60):
        xm, n = mpmath.mpf(x), a + b - 1
        exact = mpmath.fsum(
            mpmath.binomial(n, j) * xm**j * (1 - xm) ** (n - j) for j in range(a, n + 1)
        )
    with mpmath.workdps(50):
        assert close(functions._mp_betainc(*mp(a, b, x)), exact)


@pytest.mark.parametrize(
    "a, b, x",
    [
        (0.1, 0.5, 0.3),
        (0.5, 0.1, 0.99),
        (5, 0.1, 0.9),
        (0.1, 100, 0.01),
        (10, 2.5, 0.7),
        (3e4, 0.0125, 0.99998),  # small b, x near 1 (betainc_switch)
        (2900, 1000, 0.6677528530633166),  # between two switch points (betainc_switch)
    ],
)
def test_betainc_series_matches_mpmath_betainc(a, b, x):
    # Non-integer parameters, both branches of the symmetry switch.
    with mpmath.workdps(50):
        a, b, x = mp(a, b, x)
        exact = mpmath.betainc(a, b, 0, x, regularized=True)
        assert close(functions._mp_betainc(a, b, x), exact)


def test_kv_matches_upward_recurrence():
    # K_{n+1}(x) = K_{n-1}(x) + (2n/x) K_n(x) (DLMF 10.29.1), stable upwards,
    # started from K_0 and K_1. Checks the pilot point K_100(0.0631).
    with mpmath.workdps(80):
        x = mpmath.mpf(0.0631)
        prev, cur = mpmath.besselk(0, x), mpmath.besselk(1, x)
        for n in range(1, 100):
            prev, cur = cur, prev + 2 * n / x * cur
    with mpmath.workdps(50):
        value = mpmath.besselk(100, mpmath.mpf(0.0631))
        assert close(value, cur, digits=35)
        assert math.isfinite(float(value)) and float(value) > 5e305


@pytest.mark.parametrize("x", [1e-3, 1.0, 20.0, 1e3, 1e6])
def test_yv_satisfies_the_wronskian(x):
    # J_{n+1}(x) Y_n(x) - J_n(x) Y_{n+1}(x) = 2 / (pi x) (DLMF 10.5.2 with 10.6.2), for
    # every order on the grid; checks besselj against bessely across both datasets.
    with mpmath.workdps(60):
        x = mpmath.mpf(x)
        for n in (0, 1, 2, 5, 10, 20, 50, 100, 200, 500):
            w = mpmath.besselj(n + 1, x) * mpmath.bessely(n, x) - mpmath.besselj(
                n, x
            ) * mpmath.bessely(n + 1, x)
            assert close(w, 2 / (mpmath.pi * x), digits=40), n


def test_yv_matches_upward_recurrence():
    # Y_{n+1}(x) = (2n/x) Y_n(x) - Y_{n-1}(x) (DLMF 10.6.1), stable upwards for Y,
    # started from Y_0 and Y_1; reaches Y_500 at x = 20 (finite, about -1e646).
    with mpmath.workdps(80):
        x = mpmath.mpf(20)
        prev, cur = mpmath.bessely(0, x), mpmath.bessely(1, x)
        for n in range(1, 500):
            prev, cur = cur, 2 * n / x * cur - prev
    with mpmath.workdps(50):
        assert close(mpmath.bessely(500, mpmath.mpf(20)), cur, digits=35)


@pytest.mark.parametrize("nu", [0, 1, 5, 100, 500])
def test_yv_matches_hankel_expansion_at_large_x(nu):
    # Y_nu(x) ~ sqrt(2/(pi x)) (sin w sum_k (-1)^k a_2k / x^2k + cos w sum_k (-1)^k
    # a_2k+1 / x^(2k+1)), w = x - nu pi/2 - pi/4, a_k = prod_{j=1..k} (4 nu^2 - (2j-1)^2)
    # / (k! 8^k) (DLMF 10.17.4); at x = 1e6 the terms fall below 1e-45 quickly.
    with mpmath.workdps(60):
        x = mpmath.mpf(1e6)
        mu = 4 * mpmath.mpf(nu) ** 2
        p, q, term, k = mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(1), 0
        while abs(term) > mpmath.mpf(10) ** -55:
            if k % 4 == 0:
                p += term
            elif k % 4 == 1:
                q += term
            elif k % 4 == 2:
                p -= term
            else:
                q -= term
            k += 1
            term *= (mu - (2 * k - 1) ** 2) / (k * 8 * x)
        w = x - nu * mpmath.pi / 2 - mpmath.pi / 4
        expansion = mpmath.sqrt(2 / (mpmath.pi * x)) * (mpmath.sin(w) * p + mpmath.cos(w) * q)
        assert close(mpmath.bessely(nu, x), expansion, digits=40)


@pytest.mark.parametrize(
    "spec, grid, bound",
    [
        (functions.YV, functions._yv_grid, 430),
        (functions.YV_LARGE, functions._yv_large_grid, 2.2e6),
    ],
)
def test_yv_condition_numbers_match_the_tolerance_reasons(spec, grid, bound):
    # The tolerances of yv and yv_large are derived from the largest condition number
    # |x Y'/Y| at a finite point of the grid; the reasons quote these bounds.
    with mpmath.workdps(30):
        worst = max(
            functions.yv_condition(nu, x)
            for nu, x in grid()
            if abs(mpmath.bessely(nu, x)) < mpmath.mpf(2) ** 1024
        )
    assert bound * 0.97 < worst <= bound
    assert spec.tolerance >= 2 * worst * 2.0**-53


def test_grids_have_no_duplicates():
    for spec in functions.ALL:
        grid = list(spec.grid())
        assert len(grid) == len(set(grid)), spec.name


def test_betainc_series_symmetric_point_at_large_parameters():
    # I_{1/2}(a, a) = 1/2 exactly (symmetry), at the largest parameters of the
    # betainc_large grid, where the series needs the most terms.
    with mpmath.workdps(50):
        assert close(functions._mp_betainc(*mp(1e5, 1e5, 0.5)), mpmath.mpf(0.5))


@pytest.mark.parametrize("b, y", [(0.01, 1e-10), (0.02, 1e-6), (0.5, 1e-10), (0.2, 1e-1)])
def test_betainc_series_matches_closed_form_for_a_1(b, y):
    # I_x(1, b) = 1 - (1-x)^b (integrate the definition with a = 1), at x = 1 - y as in
    # betainc_switch: small b and x close to 1.
    x = 1 - y
    with mpmath.workdps(50):
        exact = 1 - (1 - mpmath.mpf(x)) ** b
        assert close(functions._mp_betainc(*mp(1, b, x)), exact)


def test_betainc_switch_grid_targets_small_b_and_the_switch_region():
    grid = functions._betainc_switch_grid()
    assert len(grid) == 186
    assert all(0 < x < 1 for _, _, x in grid)
    assert not set(grid) & set(functions._betainc_grid())
    # Small b with x close to 1, where the regression of finding 001 was.
    assert sum(b < 0.1 and x > 0.9999 for _, b, x in grid) >= 50


def test_betainc_large_grid_spans_both_tails_above_puruspe_switch():
    grid = functions._betainc_large_grid()
    assert len(grid) == 16 * 9
    assert all(a > 3000 and b > 3000 and 0 < x < 1 for a, b, x in grid)
