"""Each mpmath route used as the arbiter, checked against an independent one.

Agreement between 50 and 100 digits does not catch a wrong formula, so every
dataset's mpmath function is compared here with a different method (an exact
identity, a finite sum, or a different mpmath function) where both are valid.
"""

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


def test_grids_have_no_duplicates():
    for spec in functions.ALL:
        grid = list(spec.grid())
        assert len(grid) == len(set(grid)), spec.name
