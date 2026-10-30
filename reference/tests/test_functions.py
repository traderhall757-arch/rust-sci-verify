"""Each mpmath route used as the arbiter, checked against an independent one.

Agreement between 50 and 100 digits does not catch a wrong formula, so every
dataset's mpmath function is compared here with a different method (an exact
identity, a finite sum, or a different mpmath function) where both are valid.
"""

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

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


def test_jv_matches_power_series():
    # J_n(x) = (x/2)^n sum_k (-x^2/4)^k / (k! (n+k)!) (DLMF 10.2.2), summed at 120 digits
    # (the terms reach about 1e7 at x = 20 before they fall, so 50 digits would do).
    for n, x in [(0, 20.0), (1, 10.0), (5, 0.1), (100, 20.0), (500, 20.0), (500, 1e-10)]:
        with mpmath.workdps(120):
            xm = mpmath.mpf(x)
            q, term, k, total = -(xm**2) / 4, 1 / mpmath.factorial(n), 0, mpmath.mpf(0)
            while k < 5 or abs(term) > abs(total) * mpmath.mpf(10) ** -110:
                total += term
                k += 1
                term *= q / (k * (n + k))
            series = (xm / 2) ** n * total
        with mpmath.workdps(50):
            assert close(mpmath.besselj(n, mpmath.mpf(x)), series, digits=40), (n, x)


@pytest.mark.parametrize("nu", [0, 1, 5, 100, 500])
def test_jv_matches_hankel_expansion_at_large_x(nu):
    # J_nu(x) ~ sqrt(2/(pi x)) (cos w sum_k (-1)^k a_2k / x^2k - sin w sum_k (-1)^k
    # a_2k+1 / x^(2k+1)), with w and a_k as for Y (DLMF 10.17.3).
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
        expansion = mpmath.sqrt(2 / (mpmath.pi * x)) * (mpmath.cos(w) * p - mpmath.sin(w) * q)
        assert close(mpmath.besselj(nu, x), expansion, digits=40)


@pytest.mark.parametrize(
    "spec, grid, bound",
    [
        (functions.JV, functions._yv_grid, 200),
        (functions.JV_LARGE, functions._yv_large_grid, 2.2e6),
    ],
)
def test_jv_condition_numbers_match_the_tolerance_reasons(spec, grid, bound):
    # As for yv: the largest |x J'/J| at a normal (not underflowed) point of the grid.
    with mpmath.workdps(30):
        worst = max(
            functions.jv_condition(nu, x)
            for nu, x in grid()
            if abs(mpmath.besselj(nu, x)) >= mpmath.mpf(2) ** -1022
        )
    assert bound * 0.97 < worst <= bound
    assert spec.tolerance >= 2 * worst * 2.0**-53


def test_iv_matches_power_series():
    # I_nu(x) = (x/2)^nu sum_k (x^2/4)^k / (k! Gamma(nu+k+1)) (DLMF 10.25.2), all terms
    # positive; integer and non-integer orders from both grids.
    cases = [(0, 20.0), (1, 10.0), (5, 0.1), (100, 20.0), (500, 20.0), (500, 1e-10)]
    cases += [(0.25, 1e-5), (2.7, 5.0), (100.5, 0.1), (200.75, 20.0)]
    for nu, x in cases:
        with mpmath.workdps(80):
            xm, num = mpmath.mpf(x), mpmath.mpf(nu)
            q, term, k, total = xm**2 / 4, 1 / mpmath.gamma(num + 1), 0, mpmath.mpf(0)
            while term > total * mpmath.mpf(10) ** -70:
                total += term
                k += 1
                term *= q / (k * (num + k))
            series = (xm / 2) ** num * total
        with mpmath.workdps(50):
            assert close(mpmath.besseli(nu, mpmath.mpf(x)), series, digits=40), (nu, x)


@pytest.mark.parametrize("x", [1e-5, 0.1, 2.0, 700.0])
def test_iv_half_integer_orders_match_closed_forms(x):
    # I_{1/2}(x) = sqrt(2/(pi x)) sinh x and I_{3/2}(x) = sqrt(2/(pi x)) (cosh x - sinh x / x)
    # (DLMF 10.49.8 with 10.47.7), both on the iv_nonint grid.
    with mpmath.workdps(60):
        x = mpmath.mpf(x)
        r = mpmath.sqrt(2 / (mpmath.pi * x))
        assert close(mpmath.besseli(0.5, x), r * mpmath.sinh(x))
        assert close(mpmath.besseli(1.5, x), r * (mpmath.cosh(x) - mpmath.sinh(x) / x))


@pytest.mark.parametrize("x", [1e-3, 1.0, 20.0, 700.0])
def test_iv_satisfies_the_wronskian(x):
    # I_nu(x) K_{nu+1}(x) + I_{nu+1}(x) K_nu(x) = 1/x (DLMF 10.28.2), for the integer and
    # non-integer orders of the grids; checks besseli against besselk.
    with mpmath.workdps(60):
        x = mpmath.mpf(x)
        for nu in (0, 1, 5, 50, 500, 0.25, 2.7, 200.75):
            w = mpmath.besseli(nu, x) * mpmath.besselk(nu + 1, x)
            w += mpmath.besseli(nu + 1, x) * mpmath.besselk(nu, x)
            assert close(w, 1 / x, digits=40), nu


@pytest.mark.parametrize("nu", [0, 1, 100, 200.75])
def test_iv_matches_large_x_expansion(nu):
    # I_nu(x) ~ e^x / sqrt(2 pi x) sum_k (-1)^k a_k / x^k (DLMF 10.40.1), with a_k as for
    # J and Y. The exponentially small part this leaves out is about e^-2x = 1e-869
    # relative to the result at x = 1e3.
    with mpmath.workdps(60):
        x = mpmath.mpf(1000)
        mu = 4 * mpmath.mpf(nu) ** 2
        total, term, k = mpmath.mpf(0), mpmath.mpf(1), 0
        while abs(term) > mpmath.mpf(10) ** -55:
            total += term
            k += 1
            term *= -(mu - (2 * k - 1) ** 2) / (k * 8 * x)
        expansion = mpmath.exp(x) / mpmath.sqrt(2 * mpmath.pi * x) * total
        assert close(mpmath.besseli(nu, x), expansion, digits=40)


@pytest.mark.parametrize(
    "spec, grid, bound",
    [
        (functions.IV, functions._yv_grid, 201),
        (functions.IV_LARGE, functions._iv_large_grid, 877),
        (functions.IV_NONINT, functions._iv_nonint_grid, 728),
    ],
)
def test_iv_condition_numbers_match_the_tolerance_reasons(spec, grid, bound):
    # As for jv: the largest |x I'/I| at a normal, finite point of the grid.
    with mpmath.workdps(30):
        worst = max(
            functions.iv_condition(nu, x)
            for nu, x in grid()
            if mpmath.mpf(2) ** -1022 <= mpmath.besseli(nu, x) < mpmath.mpf(2) ** 1024
        )
    assert bound * 0.97 < worst <= bound
    assert spec.tolerance >= 2 * worst * 2.0**-53


def test_iv_large_overflows_where_the_grid_description_says():
    with mpmath.workdps(30):
        over = {(nu, x) for nu, x in functions._iv_large_grid() if mpmath.besseli(nu, x) >= 2**1024}
    expected = {(nu, 720.0) for nu in (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0)}
    expected |= {(nu, 1000.0) for nu in functions._YV_NU}
    assert over == expected


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


@pytest.mark.parametrize("x", [-38.0, -20.0, -1.5, 0.0, 0.25, 3.0, 8.0])
def test_norm_cdf_matches_ncdf_and_incomplete_gamma(x):
    # mpmath.ncdf is a separate implementation; for x < 0 also
    # Phi(x) = Gamma(1/2, x^2 / 2) / (2 sqrt(pi)) (DLMF 7.11.2 with erfc = 2 Phi(-x sqrt 2)).
    with mpmath.workdps(50):
        (x,) = mp(x)
        assert close(functions._mp_norm_cdf(x), mpmath.ncdf(x))
        if x < 0:
            route = mpmath.gammainc(mpmath.mpf(0.5), x * x / 2) / (2 * mpmath.sqrt(mpmath.pi))
            assert close(functions._mp_norm_cdf(x), route)


def test_norm_cdf_and_sf_sum_to_one_and_mirror():
    with mpmath.workdps(50):
        for x in mp(-30.0, -2.0, 0.0, 1e-300, 5.0, 1e7):
            assert close(functions._mp_norm_cdf(x) + functions._mp_norm_sf(x), 1)
            assert functions._mp_norm_sf(x) == functions._mp_norm_cdf(-x)


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-20, 0.01, 0.25, 0.5 - 1e-15, 0.75, 1 - 1e-16])
def test_norm_ppf_inverts_the_cdf(p):
    # Phi(ppf(p)) = p to far beyond double precision (relative to min(p, 1 - p)).
    with mpmath.workdps(50):
        (p,) = mp(p)
        x = functions._mp_norm_ppf(p)
        if p < 0.5:
            assert close(functions._mp_norm_cdf(x), p)
        else:
            assert close(functions._mp_norm_sf(x), 1 - p)


@pytest.mark.parametrize("p", [0.001, 0.1, 0.3, 0.5 + 1e-10, 0.9])
def test_norm_ppf_matches_inverse_error_function(p):
    # ppf(p) = sqrt(2) erfinv(2p - 1) (DLMF 7.17), a different mpmath algorithm.
    with mpmath.workdps(50):
        (p,) = mp(p)
        route = mpmath.sqrt(2) * mpmath.erfinv(2 * p - 1)
        assert close(functions._mp_norm_ppf(p), route, digits=35)
        assert functions._mp_norm_ppf(mpmath.mpf(0.5)) == 0


def test_norm_ppf_poles_are_tagged_by_rule():
    assert functions._norm_ppf_tags(0.0) == [POLE]
    assert functions._norm_ppf_tags(1.0) == [POLE]
    assert functions._norm_ppf_tags(5e-324) == []
    assert functions._norm_ppf_tags(1 - 2.0**-53) == []


@pytest.mark.parametrize("spec", [functions.NORM_CDF, functions.NORM_SF])
def test_norm_condition_numbers_match_the_tolerance_reasons(spec):
    # The largest |x phi / Phi| at a normal, nonzero point of the grid (the lower
    # tail of the cdf, mirrored for the sf).
    sign = 1 if spec is functions.NORM_CDF else -1
    with mpmath.workdps(30):
        worst = max(
            functions.norm_cdf_condition(sign * x)
            for (x,) in spec.grid()
            if spec.mpmath_fn(mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
    assert 1408 * 0.97 < worst <= 1408
    assert spec.tolerance >= 2 * worst * 2.0**-53


def test_norm_cdf_underflows_where_the_grid_description_says():
    with mpmath.workdps(30):
        cdf = functions._mp_norm_cdf
        assert cdf(mpmath.mpf(-37.51)) > mpmath.mpf(2) ** -1022 > cdf(mpmath.mpf(-37.53))
        assert cdf(mpmath.mpf(-38.48)) > mpmath.mpf(2) ** -1075 > cdf(mpmath.mpf(-38.50))
        assert functions._mp_norm_sf(mpmath.mpf(8.28)) > mpmath.mpf(2) ** -54
        assert functions._mp_norm_sf(mpmath.mpf(8.30)) < mpmath.mpf(2) ** -54


@pytest.mark.parametrize(
    ("a", "b", "x"),
    [(0.5, 0.5, 0.3), (2.0, 0.5, 0.9), (0.5, 50.0, 0.01), (25.0, 0.5, 0.97), (5000.0, 0.5, 0.9999)],
)
def test_betainc_continued_fraction_matches_series_and_mpmath(a, b, x):
    # Both branches (direct and complement) of the continued fraction against
    # the positive-term series used for betainc and against mpmath.betainc.
    with mpmath.workdps(50):
        a, b, x = mp(a, b, x)
        cf = functions._mp_betainc_cf(a, b, x)
        assert close(cf, functions._mp_betainc(a, b, x))
        assert close(cf, mpmath.betainc(a, b, 0, x, regularized=True), digits=35)


@pytest.mark.parametrize(
    "x", [-1e300, -1e5, -40.0, -3.0, -0.5, -1e-10, 0.0, 1e-300, 0.25, 7.5, 1e20]
)
def test_t_cdf_matches_closed_forms_for_df_1_2_3(x):
    # A&S 26.7.3 / 26.7.4 for df = 1 (Cauchy), 2 and 3, written for x < 0
    # without cancellation; the upper half follows from F(x) = 1 - F(-x).
    #   df = 1: F(x) = atan(-1/x) / pi
    #   df = 2: F(x) = 1 / (s (s - x)) with s = sqrt(2 + x^2)
    #   df = 3: F(x) = (2 phi - sin 2 phi) / (2 pi), phi = atan(-sqrt(3) / x)
    with mpmath.workdps(60):
        (x,) = mp(x)
        lower = -abs(x)
        if lower == 0:
            exact = {1: mpmath.mpf(1) / 2, 2: mpmath.mpf(1) / 2, 3: mpmath.mpf(1) / 2}
        else:
            s = mpmath.sqrt(2 + lower * lower)
            phi = mpmath.atan(-mpmath.sqrt(3) / lower)
            exact = {
                1: mpmath.atan(-1 / lower) / mpmath.pi,
                2: 1 / (s * (s - lower)),
                3: (2 * phi - mpmath.sin(2 * phi)) / (2 * mpmath.pi),
            }
        for df, tail in exact.items():
            if df == 3 and abs(x) > 1e5:
                continue  # 2 phi - sin 2 phi cancels beyond 60 digits
            want = tail if x <= 0 else 1 - tail
            assert close(functions._mp_t_cdf(mpmath.mpf(df), x), want), df


@pytest.mark.parametrize(
    ("df", "x"),
    [
        (0.1, -40.0),
        (0.5, 0.25),
        (2.5, -3.0),
        (7.5, -1.5),
        (50.0, 2.0),
        (1e8, -3.0),
        (1e8, -30.0),
        (1e5, -20.0),
    ],
)
def test_t_cdf_matches_hypergeometric_series(df, x):
    # F(x) = 1/2 + x Gamma((df+1)/2) / (sqrt(pi df) Gamma(df/2))
    #        * 2F1(1/2, (df+1)/2; 3/2; -x^2/df)                      (A&S 26.7.1 integrated)
    # a different mpmath algorithm (hyp2f1 with its transformations). 1/2 + ...
    # cancels to F, which is 4.9e-198 at (1e8, -30), so it runs at 300 digits.
    with mpmath.workdps(300):
        df, x = mp(df, x)
        c = mpmath.gamma((df + 1) / 2) / (mpmath.sqrt(mpmath.pi * df) * mpmath.gamma(df / 2))
        route = mpmath.mpf(1) / 2 + x * c * mpmath.hyp2f1(0.5, (df + 1) / 2, 1.5, -x * x / df)
    with mpmath.workdps(50):
        assert close(functions._mp_t_cdf(df, x), route)


@pytest.mark.parametrize(
    ("df", "x"), [(1e4, -30.0), (1e3, -40.0), (1e6, -1e3), (10.0, -1e5), (0.5, -1e300)]
)
def test_t_lower_tail_matches_beta_series_in_the_deep_tails(df, x):
    # The continued fraction against the positive-term series DLMF 8.17.8
    # (_mp_betainc, itself tested above) at h = df / (df + x^2). The series
    # needs about 1 / (1 - h) terms, so it is used only where h is not near 1.
    with mpmath.workdps(50):
        df, x = mp(df, x)
        series = functions._mp_betainc(df / 2, mpmath.mpf(1) / 2, df / (df + x * x)) / 2
        assert close(functions._mp_t_lower_tail(df, x), series)


def test_t_cdf_and_sf_sum_to_one_and_mirror():
    with mpmath.workdps(50):
        for df in mp(0.1, 3.0, 1e8):
            for x in mp(-1e300, -2.0, 0.0, 1e-10, 40.0):
                assert close(functions._mp_t_cdf(df, x) + functions._mp_t_sf(df, x), 1)
                assert functions._mp_t_sf(df, x) == functions._mp_t_cdf(df, -x)


@pytest.mark.parametrize(
    ("df", "p"),
    [
        (0.1, 1e-20),
        (0.5, 0.4999999999),
        (2.5, 5e-324),
        (10.0, 0.05),
        (1e3, 1e-100),
        (1e8, 0.25),
        (1e8, 0.999),
    ],
)
def test_t_ppf_inverts_the_cdf(df, p):
    # F(ppf(p)) = p to far beyond double precision (relative to min(p, 1 - p)).
    with mpmath.workdps(50):
        df, p = mp(df, p)
        x = functions._mp_t_ppf(df, p)
        if p < 0.5:
            assert close(functions._mp_t_cdf(df, x), p)
        else:
            assert close(functions._mp_t_sf(df, x), 1 - p)


@pytest.mark.parametrize("p", [1e-300, 1e-10, 0.25, 0.5 - 1e-10, 0.5, 0.75, 1 - 1e-15])
def test_t_ppf_matches_closed_forms_for_df_1_and_2(p):
    # df = 1: tan(pi (p - 1/2)); df = 2: (2p - 1) / sqrt(2 p (1 - p)) (A&S 26.7.1 inverted).
    # p - 1/2 and 1 - p are exact for the double p at 60 digits.
    with mpmath.workdps(60):
        (p,) = mp(p)
        cauchy = -1 / mpmath.tan(mpmath.pi * p) if p < 0.25 else mpmath.tan(mpmath.pi * (p - 0.5))
        two = (2 * p - 1) / mpmath.sqrt(2 * p * (1 - p))
        for df, want in ((1, cauchy), (2, two)):
            got = functions._mp_t_ppf(mpmath.mpf(df), p)
            assert got == want == 0 if p == 0.5 else close(got, want, digits=40)


def test_t_ppf_poles_are_tagged_by_rule():
    assert functions._t_ppf_tags(3.0, 0.0) == [POLE]
    assert functions._t_ppf_tags(3.0, 1.0) == [POLE]
    assert functions._t_ppf_tags(3.0, 5e-324) == []
    assert functions._t_ppf_tags(3.0, 1 - 2.0**-53) == []


@pytest.mark.parametrize("spec", [functions.T_CDF, functions.T_SF])
def test_t_condition_numbers_match_the_tolerance_reasons(spec):
    # The largest |x f / F| at a normal, nonzero point of the grid (the lower
    # tail of the cdf, mirrored for the sf), at df = 1e8, x = -30.
    sign = 1 if spec is functions.T_CDF else -1
    with mpmath.workdps(30):
        worst, where = max(
            (functions.t_cdf_condition(df, sign * x), (df, sign * x))
            for df, x in spec.grid()
            if spec.mpmath_fn(mpmath.mpf(df), mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
    assert 900 < worst < 902
    assert where == (1e8, -30.0)
    # 2e-13 is 2x the 1.0e-13 bound to two digits (2 * 901 * 2^-53 = 2.0006e-13).
    assert spec.tolerance >= 1.99 * worst * 2.0**-53


def test_t_ppf_log_magnitude_matches_the_tolerance_reason():
    # The largest |log|x|| at a finite, normal quantile of the committed data
    # (df = 1, p = 1e-300), against the tolerance's 2 |L| 2^-53 bound.
    data = Path(__file__).resolve().parent.parent / "data" / "t_ppf.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(abs(y))), r["in"]["df"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= abs(y) < math.inf
    ]
    worst, df, p = max(logs)
    assert (df, p) == (1.0, 1e-300)
    assert 689.6 < worst < 689.7
    assert functions.T_PPF.tolerance >= 2 * worst * 2.0**-53


@pytest.mark.parametrize(
    ("a", "x"),
    [
        (0.01, 1e-5),
        (0.01, 0.5),
        (0.5, 0.3),
        (1.5, 2.0),
        (3.0, 10.0),
        (20.0, 20.5),
        (100.0, 90.0),
        (100.0, 130.0),
        (1000.0, 500.0),
        (1000.0, 1948.6832980505137),
        (0.1, 750.0),
    ],
)
def test_gamma_cdf_and_sf_match_mpmath_gammainc(a, x):
    # Both routes (series below a + 1, continued fraction above) against
    # mpmath's own hypergeometric evaluation, where that converges.
    with mpmath.workdps(50):
        a, x = mp(a, x)
        assert close(functions._mp_gamma_cdf(a, x), mpmath.gammainc(a, 0, x, regularized=True))
        assert close(
            functions._mp_gamma_sf(a, x), mpmath.gammainc(a, x, mpmath.inf, regularized=True)
        )


@pytest.mark.parametrize(
    ("a", "x"),
    [(1, 1e-10), (2, 0.5), (5, 750.0), (1000, 970.0), (10**5, 99000.0), (10**6, 990000.0)],
)
def test_gamma_cdf_matches_finite_sum_for_integer_a(a, x):
    # Q(a, x) = e^-x sum_{k<a} x^k / k! (DLMF 8.4.10), a different method that
    # also reaches a = 1e6, where mpmath.gammainc does not converge and SciPy
    # is 2.1e-10 off (P = 5.4e-24 at x = 990000). P = 1 - Q cancels to P, so
    # the sum runs with 60 digits more than P's magnitude needs.
    xm = mpmath.mpf(x)
    with mpmath.workdps(110):
        term = total = mpmath.mpf(1)
        for k in range(1, a):
            term *= xm / k
            total += term
        q = mpmath.exp(-xm) * total
        p = 1 - q
    with mpmath.workdps(50):
        assert close(functions._mp_gamma_sf(mpmath.mpf(a), xm), q)
        assert close(functions._mp_gamma_cdf(mpmath.mpf(a), xm), p)


@pytest.mark.parametrize("x", [1e-300, 1e-10, 0.3, 1.0, 2.0, 30.0, 700.0, 1e5])
def test_gamma_cdf_matches_closed_forms_for_a_one_half_and_one(x):
    # a = 1/2: P = erf(sqrt x), Q = erfc(sqrt x); a = 1: P = 1 - e^-x, Q = e^-x.
    with mpmath.workdps(50):
        (x,) = mp(x)
        half = mpmath.mpf(1) / 2
        assert close(functions._mp_gamma_cdf(half, x), mpmath.erf(mpmath.sqrt(x)))
        assert close(functions._mp_gamma_sf(half, x), mpmath.erfc(mpmath.sqrt(x)))
        assert close(functions._mp_gamma_cdf(mpmath.mpf(1), x), -mpmath.expm1(-x))
        assert close(functions._mp_gamma_sf(mpmath.mpf(1), x), mpmath.exp(-x))


def test_gamma_cdf_and_sf_sum_to_one_and_cover_the_ends():
    with mpmath.workdps(50):
        for a in mp(0.01, 3.0, 1e6):
            assert functions._mp_gamma_cdf(a, mpmath.mpf(0)) == 0
            assert functions._mp_gamma_sf(a, mpmath.mpf(0)) == 1
            for x in mp(1e-300, 0.5, 1e6 - 1, 1e6 + 1, 1e300):
                total = functions._mp_gamma_cdf(a, x) + functions._mp_gamma_sf(a, x)
                assert close(total, 1)


@pytest.mark.parametrize(
    ("a", "p"),
    [
        (0.01, 1e-3),
        (0.1, 5e-324),
        (0.5, 0.4999999999),
        (3.0, 0.05),
        (20.0, 1 - 1e-15),
        (1e3, 1e-100),
        (1e3, 5e-324),
        (1e5, 5e-324),
        (1e6, 1e-10),
        (1e6, 0.25),
        (1e6, 0.999),
    ],
)
def test_gamma_ppf_inverts_the_cdf(a, p):
    # P(ppf(p)) = p, or Q(ppf(p)) = 1 - p in the upper half, far beyond double
    # precision; includes the points where SciPy is off (p = 5e-324 for
    # a >= 100, and a = 1e6, p = 1e-10).
    with mpmath.workdps(50):
        a, p = mp(a, p)
        x = functions._mp_gamma_ppf(a, p)
        if p <= 0.5:
            assert close(functions._mp_gamma_cdf(a, x), p)
        else:
            assert close(functions._mp_gamma_sf(a, x), 1 - p)


@pytest.mark.parametrize("p", [1e-300, 1e-10, 0.25, 0.5, 0.75, 1 - 1e-15])
def test_gamma_ppf_matches_closed_forms_for_a_one_half_and_one(p):
    # a = 1: x = -log(1 - p); a = 1/2: x = erfinv(p)^2 (p = erf(sqrt x)); in the
    # upper half x = erfcinv(1 - p)^2, written with erfinv to stay accurate.
    # 1 - p is exact for the double p at 60 digits.
    with mpmath.workdps(60):
        (p,) = mp(p)
        assert close(functions._mp_gamma_ppf(mpmath.mpf(1), p), -mpmath.log1p(-p))
        half = functions._mp_gamma_ppf(mpmath.mpf(1) / 2, p)
        assert close(half, mpmath.erfinv(p) ** 2, digits=35)
    assert functions._mp_gamma_ppf(mpmath.mpf(3), mpmath.mpf(0)) == 0


def test_gamma_ppf_poles_are_tagged_by_rule():
    assert functions._gamma_ppf_tags(3.0, 1.0) == [POLE]
    assert functions._gamma_ppf_tags(3.0, 0.0) == []
    assert functions._gamma_ppf_tags(3.0, 5e-324) == []
    assert functions._gamma_ppf_tags(3.0, 1 - 2.0**-53) == []


@pytest.mark.parametrize(
    ("spec", "condition", "where"),
    [
        (functions.GAMMA_CDF, functions.gamma_cdf_condition, (1e6, 1e6 - 30 * 1e3)),
        (functions.GAMMA_SF, functions.gamma_sf_condition, (1e6, 1e6 + 30 * 1e3)),
    ],
)
def test_gamma_condition_numbers_match_the_tolerance_reasons(spec, condition, where):
    # The largest condition number at a normal, nonzero point of the grid is
    # 3.0e4 at a = 1e6, x = a -+ 30 sqrt(a): what a method that rounds x / a
    # loses there (3.3e-12) would not fit the tolerance, as the reason says.
    with mpmath.workdps(30):
        worst, at = max(
            (condition(a, x), (a, x))
            for a, x in spec.grid()
            if x > 0 and spec.mpmath_fn(mpmath.mpf(a), mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
    assert 3.0e4 < worst < 3.01e4
    assert at == where
    assert 3.3e-12 < worst * 2.0**-53 < 3.4e-12
    assert worst * 2.0**-53 > spec.tolerance
    # A careful method's k^2 2^-53 at k = 30 fits with room to spare.
    assert 30**2 * 2.0**-53 < spec.tolerance / 9


def test_gamma_ppf_log_magnitude_matches_the_tolerance_reason():
    # The largest |log x| at a nonzero, normal quantile of the committed data
    # (a = 0.01, p = 1e-3, x = 5.7e-301), against the tolerance's 2 |L| 2^-53 bound.
    data = Path(__file__).resolve().parent.parent / "data" / "gamma_ppf.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y)), r["in"]["a"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    worst, a, p = max(logs)
    assert (a, p) == (0.01, 1e-3)
    assert 691.3 < worst < 691.4
    assert functions.GAMMA_PPF.tolerance >= 2 * worst * 2.0**-53


@pytest.mark.parametrize(
    ("a", "x"),
    [(0.01, 1e-3), (0.5, 0.5), (1.0, 1 / 750), (5.0, 1 / 750), (50.0, 1e5), (1e3, 1 / 1948.6),
     (1e4, 1e-4), (1e4, 1 / 13000), (3.0, 1e300)],
)  # fmt: skip
def test_invgamma_cdf_and_sf_match_mpmath_gammainc_at_the_reciprocal(a, x):
    # F(x) = Q(a, 1/x) and 1 - F(x) = P(a, 1/x), by mpmath.gammainc (a different
    # method from the series / continued fraction behind the datasets) with 1/x
    # formed at 60 extra bits; includes the subnormal results SciPy flushes to 0
    # (a = 5, x = 1/750 and a = 50, x = 1e5).
    with mpmath.workdps(50):
        a, x = mp(a, x)
        with mpmath.extraprec(60):
            t = 1 / x
        q = mpmath.gammainc(a, t, mpmath.inf, regularized=True)
        p = mpmath.gammainc(a, 0, t, regularized=True)
        assert close(functions._mp_invgamma_cdf(a, x), q)
        assert close(functions._mp_invgamma_sf(a, x), p)
        assert close(functions._mp_invgamma_cdf(a, x) + functions._mp_invgamma_sf(a, x), 1)


@pytest.mark.parametrize("x", [1e-3, 1 / 750, 0.5, 1.0, 3.0, 1e5, 1e300])
def test_invgamma_cdf_and_sf_match_closed_forms_for_a_one_half_and_one(x):
    # a = 1: F = e^(-1/x), 1 - F = -expm1(-1/x); a = 1/2: F = erfc(1/sqrt(x)),
    # 1 - F = erf(1/sqrt(x)).
    with mpmath.workdps(50):
        (x,) = mp(x)
        one, half = mpmath.mpf(1), mpmath.mpf(1) / 2
        assert close(functions._mp_invgamma_cdf(one, x), mpmath.exp(-1 / x))
        assert close(functions._mp_invgamma_sf(one, x), -mpmath.expm1(-1 / x))
        assert close(functions._mp_invgamma_cdf(half, x), mpmath.erfc(1 / mpmath.sqrt(x)))
        assert close(functions._mp_invgamma_sf(half, x), mpmath.erf(1 / mpmath.sqrt(x)))


def test_invgamma_cdf_and_sf_cover_the_ends():
    with mpmath.workdps(50):
        for a in mp(0.01, 3.0, 1e4):
            assert functions._mp_invgamma_cdf(a, mpmath.mpf(0)) == 0
            assert functions._mp_invgamma_sf(a, mpmath.mpf(0)) == 1
            # Where 1/x overflows a double the cdf still rounds to 0 and the sf to 1
            # (only the rounded double matters there: e^(-1/x) is about 10^(-4e309)).
            assert functions._mp_invgamma_cdf(a, mpmath.mpf(5e-324)) < mpmath.mpf(2) ** -1075
            assert 1 - functions._mp_invgamma_sf(a, mpmath.mpf(1e-310)) < mpmath.mpf(2) ** -54


@pytest.mark.parametrize(
    ("a", "p"),
    [(0.01, 5e-324), (0.1, 1 - 1e-15), (0.5, 5e-324), (0.5, 1 - 1e-15), (1.0, 1e-300),
     (3.0, 0.5), (5.0, 5e-324), (1e3, 5e-324), (1e4, 5e-324), (1e4, 0.5000000001),
     (0.01, 0.999)],
)  # fmt: skip
def test_invgamma_ppf_inverts_the_cdf(a, p):
    # F(ppf(p)) = p, or 1 - F(ppf(p)) = 1 - p in the upper half (1 - p is exact
    # at 50 digits), with F by mpmath.gammainc at 1/x; includes every kind of
    # point where SciPy or statrs is off (p = 5e-324, 1 - 1e-15 for small a).
    with mpmath.workdps(50):
        a, p = mp(a, p)
        x = functions._mp_invgamma_ppf(a, p)
        t = 1 / x
        if p <= 0.5:
            assert close(mpmath.gammainc(a, t, mpmath.inf, regularized=True), p)
        else:
            assert close(mpmath.gammainc(a, 0, t, regularized=True), 1 - p)


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-10, 0.25, 0.5, 0.75, 1 - 1e-15])
def test_invgamma_ppf_matches_the_closed_form_for_a_one(p):
    # a = 1: F(x) = e^(-1/x), so x = -1 / log(p); log(p) is exact enough for the
    # double p at 60 digits, also next to 1.
    with mpmath.workdps(60):
        (p,) = mp(p)
        assert close(functions._mp_invgamma_ppf(mpmath.mpf(1), p), -1 / mpmath.log(p))
    assert functions._mp_invgamma_ppf(mpmath.mpf(3), mpmath.mpf(0)) == 0


@pytest.mark.parametrize(
    ("spec", "condition"),
    [
        (functions.INVGAMMA_CDF, functions.invgamma_cdf_condition),
        (functions.INVGAMMA_SF, functions.invgamma_sf_condition),
    ],
)
def test_invgamma_condition_numbers_match_the_tolerance_reason(spec, condition):
    # The largest condition number at a normal, nonzero point is 3.0e3 at
    # a = 1e4, near 1/x = a -+ 30 sqrt(a): rounding 1/x costs at most 3.4e-13
    # there, inside the tolerance, while gamma's a = 1e5 and 1e6 would cost
    # 1.1e-12 and 3.3e-12, beyond it.
    with mpmath.workdps(30):
        worst, a = max(
            (condition(a, x), a)
            for a, x in spec.grid()
            if x > 0 and spec.mpmath_fn(mpmath.mpf(a), mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
        assert a == 1e4
        assert 3.0e3 < worst < 3.01e3
        assert worst * 2.0**-53 < 3.4e-13 < spec.tolerance
        for big, low, high in ((1e5, 1.0e-12, 1.1e-12), (1e6, 3.2e-12, 3.4e-12)):
            k30 = condition(big, 1 / (big + 30 * math.sqrt(big)))
            k30 = max(k30, condition(big, 1 / (big - 30 * math.sqrt(big))))
            assert low < k30 * 2.0**-53 < high
            assert k30 * 2.0**-53 > spec.tolerance


def test_invgamma_ppf_log_magnitude_matches_the_tolerance_reason():
    # The largest |log x| at a finite, normal quantile of the committed data
    # (a = 0.01, p = 0.999, x = 1.8e300), against the tolerance's 2 |L| 2^-53 bound.
    data = Path(__file__).resolve().parent.parent / "data" / "invgamma_ppf.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y)), r["in"]["a"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    worst, a, p = max(logs)
    assert (a, p) == (0.01, 0.999)
    assert 691.3 < worst < 691.4
    assert functions.INVGAMMA_PPF.tolerance >= 2 * worst * 2.0**-53


@pytest.mark.parametrize("x", [5e-324, 1e-170, 4e-162, 1e-160, 1e-3, 0.5, 1.0, 3.0, 38.0, 1e4])
def test_chi_cdf_and_sf_match_closed_forms_for_k_one_two_three(x):
    # k = 1: F = erf(x / sqrt 2); k = 2: 1 - F = e^(-x^2/2); k = 3:
    # 1 - F = erfc(x / sqrt 2) + sqrt(2 / pi) x e^(-x^2/2) (DLMF 8.4.6 and 8.4.10),
    # computed from x itself, so also where x * x is subnormal or 0 in a double
    # (mpmath.erfc overflows internally above about 1e154: see the ends test).
    with mpmath.workdps(50):
        (x,) = mp(x)
        one, two, three = mp(1, 2, 3)
        s2 = mpmath.sqrt(2)
        e = mpmath.exp(-x * x / 2)
        assert close(functions._mp_chi_cdf(one, x), mpmath.erf(x / s2))
        assert close(functions._mp_chi_sf(one, x), mpmath.erfc(x / s2))
        assert close(functions._mp_chi_cdf(two, x), -mpmath.expm1(-x * x / 2))
        assert close(functions._mp_chi_sf(two, x), e)
        sf3 = mpmath.erfc(x / s2) + mpmath.sqrt(2 / mpmath.pi) * x * e
        assert close(functions._mp_chi_sf(three, x), sf3)
        assert close(functions._mp_chi_cdf(three, x) + sf3, 1)


@pytest.mark.parametrize(
    ("k", "x"),
    [(1.0, 1e-160), (5.0, 2.0), (100.0, 0.00447213595499958), (161.0, 44.721359549995796),
     (1e3, 20.0), (1e3, 37.416573867739416)],
)  # fmt: skip
def test_chi_cdf_and_sf_match_mpmath_gammainc_at_half_the_square(k, x):
    # F(x) = P(k/2, x^2/2) by mpmath.gammainc (a different method from the series
    # and continued fraction behind the datasets), with x^2/2 formed exactly.
    with mpmath.workdps(50):
        k, x = mp(k, x)
        t = x * x / 2  # exact: x has 53 bits and the precision is 169
        assert close(functions._mp_chi_cdf(k, x), mpmath.gammainc(k / 2, 0, t, regularized=True))
        assert close(
            functions._mp_chi_sf(k, x), mpmath.gammainc(k / 2, t, mpmath.inf, regularized=True)
        )


@pytest.mark.parametrize("x", [75.87726479572596, 119.34253511267173])
def test_chi_cdf_and_sf_match_the_finite_sum_at_k_ten_thousand(x):
    # Where SciPy is 1.3e-12 (cdf) and 5.1e-12 (sf) off: Q(a, t) = e^-t
    # sum_{j<a} t^j / j! (DLMF 8.4.10) for the integer a = k/2 = 5000, with
    # 60 digits more than P = 1 - Q's magnitude needs.
    xm = mpmath.mpf(x)
    with mpmath.workdps(400):
        t = xm * xm / 2
        term = total = mpmath.mpf(1)
        for j in range(1, 5000):
            term *= t / j
            total += term
        q = mpmath.exp(-t) * total
        p = 1 - q
    with mpmath.workdps(50):
        assert close(functions._mp_chi_sf(mpmath.mpf(1e4), xm), q)
        assert close(functions._mp_chi_cdf(mpmath.mpf(1e4), xm), p)


def test_chi_cdf_and_sf_cover_the_ends():
    with mpmath.workdps(50):
        for k in mp(1, 161, 1e4):
            assert functions._mp_chi_cdf(k, mpmath.mpf(0)) == 0
            assert functions._mp_chi_sf(k, mpmath.mpf(0)) == 1
            # Where x * x overflows a double, the cdf still rounds to 1 and the sf to 0.
            for x in mp(1e155, 1e300):
                assert 1 - functions._mp_chi_cdf(k, x) < mpmath.mpf(2) ** -54
                assert functions._mp_chi_sf(k, x) < mpmath.mpf(2) ** -1075


@pytest.mark.parametrize(
    ("k", "p"),
    [(1.0, 5e-324), (1.0, 1e-300), (2.0, 5e-324), (10.0, 5e-324), (300.0, 5e-324),
     (1e4, 5e-324), (3.0, 0.5), (161.0, 0.999), (1.0, 1 - 1e-15), (1e4, 1 - 1e-15)],
)  # fmt: skip
def test_chi_ppf_inverts_the_cdf(k, p):
    # F(ppf(p)) = p, or 1 - F(ppf(p)) = 1 - p in the upper half (1 - p exact),
    # with F by mpmath.gammainc at x^2/2; includes every kind of point where
    # SciPy or statrs is off.
    with mpmath.workdps(50):
        k, p = mp(k, p)
        x = functions._mp_chi_ppf(k, p)
        t = x * x / 2
        if p <= 0.5:
            assert close(mpmath.gammainc(k / 2, 0, t, regularized=True), p, 35)
        else:
            assert close(mpmath.gammainc(k / 2, t, mpmath.inf, regularized=True), 1 - p, 35)


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-158, 1e-10, 0.25, 0.5, 0.75, 1 - 1e-15])
def test_chi_ppf_matches_closed_forms_for_k_one_and_two(p):
    # k = 1: x = sqrt(2) erfinv(p); k = 2: x = sqrt(-2 log(1 - p)).
    with mpmath.workdps(60):
        (p,) = mp(p)
        one, two = mp(1, 2)
        assert close(functions._mp_chi_ppf(one, p), mpmath.sqrt(2) * mpmath.erfinv(p))
        assert close(functions._mp_chi_ppf(two, p), mpmath.sqrt(-2 * mpmath.log1p(-p)))
    assert functions._mp_chi_ppf(mpmath.mpf(3), mpmath.mpf(0)) == 0


@pytest.mark.parametrize(
    ("spec", "condition", "gamma_condition"),
    [
        (functions.CHI_CDF, functions.chi_cdf_condition, functions.gamma_cdf_condition),
        (functions.CHI_SF, functions.chi_sf_condition, functions.gamma_sf_condition),
    ],
)
def test_chi_condition_numbers_match_the_tolerance_reason(spec, condition, gamma_condition):
    # Rounding x * x costs gamma's condition number at t = x^2/2 (half the chi
    # condition number) times 2^-53: at most 2.1e3 at k = 1e4, near
    # t = a -+ 30 sqrt(a), so 2.4e-13, inside the tolerance; at k = 2e6
    # (gamma's a = 1e6) it would cost 3.3e-12, beyond it.
    with mpmath.workdps(30):
        assert close(condition(5.0, 2.0), 2 * gamma_condition(2.5, 2.0), 25)
        worst, k = max(
            (condition(k, x) / 2, k)
            for k, x in spec.grid()
            if 0 < x < 1e154
            and spec.mpmath_fn(mpmath.mpf(k), mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
        assert k == 1e4
        assert 2.12e3 < worst < 2.13e3
        assert worst * 2.0**-53 < 2.4e-13 < spec.tolerance
        a = 1e6
        k30 = max(
            gamma_condition(a, a + 30 * math.sqrt(a)), gamma_condition(a, a - 30 * math.sqrt(a))
        )
        assert 3.3e-12 < k30 * 2.0**-53 < 3.4e-12


def test_chi_ppf_log_magnitude_matches_the_tolerance_reason():
    # The largest |log(x / sqrt 2)| at a nonzero, normal quantile of the committed
    # data (k = 1, p = 1e-300, x = 1.25e-300), against the 2 |L / 2| 2^-53 bound.
    data = Path(__file__).resolve().parent.parent / "data" / "chi_ppf.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y / math.sqrt(2))), r["in"]["k"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    worst, k, p = max(logs)
    assert (k, p) == (1.0, 1e-300)
    assert 690.8 < worst < 691.0
    assert functions.CHI_PPF.tolerance >= 2 * worst * 2.0**-53


@pytest.mark.parametrize("x", [5e-324, 1e-300, 1e-16, 1e-5, 0.5, 1.0, 37.5, 708.5, 745.0])
def test_expon_cdf_and_sf_match_the_gamma_routes_and_the_series(x):
    # The exponential distribution is the gamma distribution with shape 1:
    # F = P(1, x) (the gamma_cdf series / continued fraction, and
    # mpmath.gammainc), Q = Q(1, x). For small x also the series
    # x - x^2/2 + x^3/6 - ..., summed directly.
    with mpmath.workdps(50):
        (x,) = mp(x)
        one = mpmath.mpf(1)
        assert close(functions._mp_expon_cdf(x), functions._mp_gamma_cdf(one, x))
        assert close(functions._mp_expon_sf(x), functions._mp_gamma_sf(one, x))
        route = mpmath.gammainc(one, 0, x, regularized=True)
        assert close(functions._mp_expon_cdf(x), route)
        assert close(functions._mp_expon_cdf(x) + functions._mp_expon_sf(x), 1)
        if x < 1e-3:
            series = mpmath.fsum((-1) ** (k + 1) * x**k / mpmath.factorial(k) for k in range(1, 40))
            assert close(functions._mp_expon_cdf(x), series)


def test_expon_cdf_and_sf_below_the_support():
    with mpmath.workdps(50):
        for x in mp(-1e300, -1.0, -0.25, 0.0):
            assert functions._mp_expon_cdf(x) == 0
            assert functions._mp_expon_sf(x) == 1


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-16, 0.25, 0.5, 0.75, 1 - 1e-10, 1 - 2.0**-53])
def test_expon_ppf_inverts_the_cdf_and_matches_the_log(p):
    # F(ppf(p)) = p (relative to min(p, 1 - p)), and ppf(p) = -log(1 - p) with
    # 1 - p exact (1200 bits hold 1 - 2^-1074).
    with mpmath.workdps(50):
        (p,) = mp(p)
        x = functions._mp_expon_ppf(p)
        if p < 0.5:
            assert close(functions._mp_expon_cdf(x), p)
        else:
            assert close(functions._mp_expon_sf(x), 1 - p)
        with mpmath.workprec(1200):
            route = -mpmath.log(1 - p)
        assert close(x, route)
        assert functions._mp_expon_ppf(mpmath.mpf(0)) == 0


def test_expon_ppf_poles_are_tagged_by_rule():
    assert functions._expon_ppf_tags(1.0) == [POLE]
    assert functions._expon_ppf_tags(0.0) == []
    assert functions._expon_ppf_tags(1 - 2.0**-53) == []


def test_expon_condition_numbers_match_the_tolerance_reasons():
    # cdf: at most 1 (x / (e^x - 1)); sf: x, 708.4 at the last normal result.
    with mpmath.workdps(30):
        grid = [x for (x,) in functions.EXPON_CDF.grid() if x > 0]
        assert max(functions.expon_cdf_condition(x) for x in grid) <= 1
        last_normal = max(
            x for x in grid if functions._mp_expon_sf(mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
        assert last_normal == 708.25
        assert (
            functions.expon_sf_condition(708.4) * 2.0**-53 < 7.9e-14 < functions.EXPON_SF.tolerance
        )


def test_expon_sf_underflows_where_the_grid_description_says():
    with mpmath.workdps(30):
        sf = functions._mp_expon_sf
        assert sf(mpmath.mpf(708.39)) > mpmath.mpf(2) ** -1022 > sf(mpmath.mpf(708.41))
        assert sf(mpmath.mpf(745.13)) > mpmath.mpf(2) ** -1075 > sf(mpmath.mpf(745.14))
        assert functions._mp_expon_cdf(mpmath.mpf(37.4)) < 1 - mpmath.mpf(2) ** -54
        assert functions._mp_expon_cdf(mpmath.mpf(37.5)) > 1 - mpmath.mpf(2) ** -54


@pytest.mark.parametrize("x", [5e-324, 1e-300, 1e-16, 1e-5, 0.5, 1.0, 37.5, 708.5, 745.0])
def test_weibull_with_shape_1_is_the_exponential(x):
    # Shape 1, scale 1 is the exponential distribution with rate 1.
    with mpmath.workdps(50):
        (x,) = mp(x)
        one = mpmath.mpf(1)
        assert close(functions._mp_weibull_cdf(one, one, x), functions._mp_expon_cdf(x))
        assert close(functions._mp_weibull_sf(one, one, x), functions._mp_expon_sf(x))


@pytest.mark.parametrize(
    ("k", "scale", "x"),
    [(0.1, 2.0**-500, 2.0**-510), (0.5, 1.0, 0.3), (1.5, 1024.0, 900.0), (2.1, 2.0**500, 2.0**501),
     (3.5, 2.0**-10, 2.0**-11), (10.0, 1.0, 1.2), (100.0, 1024.0, 1030.0)],
)  # fmt: skip
def test_weibull_cdf_and_sf_match_the_integral_of_the_density(k, scale, x):
    # F(x) = integral of f(s) = (k / lambda) (s/lambda)^(k-1) exp(-(s/lambda)^k)
    # from 0 to x, by quadrature in w with s/lambda = w^m, m = ceil(1/k), which
    # removes the singularity at 0 for k < 1: the integrand is
    # m k w^(m k - 1) exp(-w^(m k)).
    with mpmath.workdps(50):
        k, scale, x = mp(k, scale, x)
        m = math.ceil(1 / k)
        end = (x / scale) ** (mpmath.mpf(1) / m)
        integral = mpmath.quad(
            lambda w: m * k * w ** (m * k - 1) * mpmath.exp(-(w ** (m * k))), [0, end / 2, end]
        )
        cdf = functions._mp_weibull_cdf(k, scale, x)
        assert close(cdf, integral)
        assert close(cdf + functions._mp_weibull_sf(k, scale, x), 1)


def test_weibull_cdf_and_sf_below_the_support():
    with mpmath.workdps(50):
        for k in mp(0.1, 2.0, 100.0):
            for x in mp(-1e300, -1.0, 0.0):
                assert functions._mp_weibull_cdf(k, mpmath.mpf(2), x) == 0
                assert functions._mp_weibull_sf(k, mpmath.mpf(2), x) == 1


def test_weibull_grid_has_exact_standardised_x():
    # Every positive x is u * lambda with u and x normal doubles, so x / lambda
    # is exact (the tolerance reasons rely on it).
    grid = functions.WEIBULL_CDF.grid()
    assert len(grid) == 1139
    for k, scale, x in grid:
        assert math.frexp(scale)[0] == 0.5  # a power of two
        if x > 0:
            u = x / scale
            assert 2.0**-1022 <= u < math.inf and u * scale == x


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-16, 0.25, 0.5, 0.75, 1 - 1e-10, 1 - 2.0**-53])
@pytest.mark.parametrize(("k", "scale"), [(0.1, 2.0**-500), (2.1, 1.0), (100.0, 2.0**500)])
def test_weibull_ppf_inverts_the_cdf_and_matches_the_log(k, scale, p):
    # F(ppf(p)) = p (relative to min(p, 1 - p)), and ppf(p) = lambda
    # (-log(1 - p))^(1/k) with 1 - p exact (1200 bits hold 1 - 2^-1074).
    with mpmath.workdps(50):
        k, scale, p = mp(k, scale, p)
        x = functions._mp_weibull_ppf(k, scale, p)
        if p < 0.5:
            assert close(functions._mp_weibull_cdf(k, scale, x), p)
        else:
            assert close(functions._mp_weibull_sf(k, scale, x), 1 - p)
        with mpmath.workprec(1200):
            route = scale * mpmath.exp(mpmath.log(-mpmath.log(1 - p)) / k)
        assert close(x, route)
        assert functions._mp_weibull_ppf(k, scale, mpmath.mpf(0)) == 0


def test_weibull_ppf_poles_are_tagged_by_rule():
    assert functions._weibull_ppf_tags(2.0, 1.0, 1.0) == [POLE]
    assert functions._weibull_ppf_tags(2.0, 1.0, 0.0) == []
    assert functions._weibull_ppf_tags(2.0, 1.0, 1 - 2.0**-53) == []


def test_weibull_condition_numbers_match_the_tolerance_reasons():
    # cdf: at most 1 in t; sf: t, below 708.4 at every normal result.
    with mpmath.workdps(30):
        grid = [(k, s, x) for k, s, x in functions.WEIBULL_CDF.grid() if x > 0]
        assert max(functions.weibull_cdf_condition_in_t(*p) for p in grid) <= 1
        normal = [p for p in grid if functions._mp_weibull_sf(*mp(*p)) >= mpmath.mpf(2) ** -1022]
        worst = max(functions.weibull_sf_condition_in_t(*p) for p in normal)
        assert 700 < worst < 708.4
        assert 708.4 * 2.0**-53 < 7.9e-14
        assert 708.4 * 1.04 * 2.0**-53 < 8.2e-14 < functions.WEIBULL_SF.tolerance


def test_weibull_ppf_rounding_of_one_over_k_matches_the_tolerance_reason():
    # Rounding 1/k costs |log(x/lambda)| 2^-53; the largest at a normal result.
    path = Path(__file__).resolve().parent.parent / "data" / "weibull_ppf.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y / r["in"]["lambda"][0])), r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    worst, p = max(logs)
    assert p == 5e-324
    assert 744.4 < worst < 744.5
    assert worst * 2.0**-53 < 8.3e-14 < functions.WEIBULL_PPF.tolerance


@pytest.mark.parametrize(("mu", "sigma"), [(0.0, 0.25), (-3.0, 0.5), (0.0, 1.0), (2.0, 4.0)])
def test_lognorm_cdf_and_sf_match_mpmath_ncdf(mu, sigma):
    # Checks the wiring (ln x, mu, sigma, the sign of w) at every grid point.
    # mpmath.ncdf is erfc(-w / sqrt 2) / 2 for w < 0, the same formula as
    # _mp_norm_cdf, and (1 + erf(w / sqrt 2)) / 2 for w >= 0; the independent
    # check of the tails is the quadrature test below.
    with mpmath.workdps(50):
        for m, s, x in functions.LOGNORM_CDF.grid():
            if (m, s) != (mu, sigma) or x <= 0:
                continue
            m, s, x = mp(m, s, x)
            w = (mpmath.log(x) - m) / s
            assert close(functions._mp_lognorm_cdf(m, s, x), mpmath.ncdf(w))
            assert close(functions._mp_lognorm_sf(m, s, x), mpmath.ncdf(-w))


@pytest.mark.parametrize(
    ("mu", "sigma", "w"),
    [(0.0, 0.25, -1.0), (-3.0, 0.5, -37.5), (-3.0, 0.5, -8.0), (0.0, 1.0, 0.0), (0.0, 1.0, 1.5),
     (0.0, 1.0, -30.0), (2.0, 4.0, -12.0), (2.0, 4.0, 37.5), (0.0, 0.25, 20.0)],
)  # fmt: skip
def test_lognorm_cdf_and_sf_match_the_integral_of_the_density(mu, sigma, w):
    # F(x) = integral of f(s) = phi((ln s - mu) / sigma) / (s sigma) from 0 to
    # x; with s = exp(mu + sigma (w - u)), w = (ln x - mu) / sigma, the
    # integrand is phi(w - u) for u from 0 to infinity, and 1 - F(x) is the
    # integral of phi(w + u). Both tails are integrated directly, out to the
    # last normal result (|w| = 37.5), so no erfc is involved. The breakpoints
    # follow the e^(-|w| u) decay; quadrature gives 47 or more digits here.
    with mpmath.workdps(50):
        x = float(mpmath.exp(mu + sigma * mpmath.mpf(w)))
        mu, sigma, x = mp(mu, sigma, x)
        w = (mpmath.log(x) - mu) / sigma
        h = 1 / (1 + abs(w))
        points = [0, h, 4 * h, 16 * h, 64 * h, mpmath.inf]

        def tail(sign):
            # phi(w + sign u) = phi(w) exp(-sign w u - u^2 / 2); phi(w) is taken
            # out so the integrand is at most 1.
            phi = mpmath.exp(-(w**2) / 2) / mpmath.sqrt(2 * mpmath.pi)
            return phi * mpmath.quad(lambda u: mpmath.exp(-sign * w * u - u**2 / 2), points)

        cdf = functions._mp_lognorm_cdf(mu, sigma, x)
        sf = functions._mp_lognorm_sf(mu, sigma, x)
        # The smaller of F and 1 - F is integrated (the larger one's integrand
        # peaks far from the breakpoints); the other follows from F + (1 - F) = 1.
        if w <= 0:
            assert close(cdf, tail(-1))
        if w >= 0:
            assert close(sf, tail(1))
        assert close(cdf + sf, 1)


def test_lognorm_cdf_and_sf_below_the_support():
    with mpmath.workdps(50):
        for x in mp(-1e300, -1.0, 0.0):
            assert functions._mp_lognorm_cdf(mpmath.mpf(0), mpmath.mpf(1), x) == 0
            assert functions._mp_lognorm_sf(mpmath.mpf(0), mpmath.mpf(1), x) == 1


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-16, 0.25, 0.5, 0.75, 1 - 1e-10, 1 - 2.0**-53])
@pytest.mark.parametrize(("mu", "sigma"), [(0.0, 0.25), (-3.0, 0.5), (2.0, 4.0)])
def test_lognorm_ppf_inverts_the_cdf_and_matches_erfinv(mu, sigma, p):
    # F(ppf(p)) = p (relative to min(p, 1 - p)), and ppf(p) = exp(mu + sigma
    # sqrt(2) erfinv(2p - 1)) with 2p - 1 exact at 1200 bits (for p >= 1e-300).
    with mpmath.workdps(50):
        mu, sigma, p = mp(mu, sigma, p)
        x = functions._mp_lognorm_ppf(mu, sigma, p)
        if p < 0.5:
            assert close(functions._mp_lognorm_cdf(mu, sigma, x), p)
        else:
            assert close(functions._mp_lognorm_sf(mu, sigma, x), 1 - p)
        if p >= 1e-300:
            with mpmath.workprec(1200):
                route = mpmath.exp(mu + sigma * mpmath.sqrt(2) * mpmath.erfinv(2 * p - 1))
            assert close(x, route)
        assert functions._mp_lognorm_ppf(mu, sigma, mpmath.mpf(0)) == 0


def test_lognorm_ppf_poles_are_tagged_by_rule():
    assert functions._lognorm_ppf_tags(0.0, 1.0, 1.0) == [POLE]
    assert functions._lognorm_ppf_tags(0.0, 1.0, 0.0) == []
    assert functions._lognorm_ppf_tags(0.0, 1.0, 1 - 2.0**-53) == []


def test_lognorm_rounding_costs_match_the_tolerance_reasons():
    # The largest cost of rounding ln x and w at a normal result, in 2^-53.
    with mpmath.workdps(30):
        grid = [p for p in functions.LOGNORM_CDF.grid() if p[2] > 0]
        tiny = mpmath.mpf(2) ** -1022
        for upper, fn, bound, at in (
            (False, functions._mp_lognorm_cdf, 3.4e-13, (-3.0, 0.5, -37.5)),
            (True, functions._mp_lognorm_sf, 3.2e-13, (2.0, 4.0, 37.5)),
        ):
            normal = [p for p in grid if fn(*mp(*p)) >= tiny]
            cost, point = max((functions.lognorm_rounding_cost(*p, upper), p) for p in normal)
            assert 0.9 * bound < cost * 2.0**-53 < bound < functions.LOGNORM_CDF.tolerance / 2
            mu, sigma, x = mp(*point)
            assert point[:2] == at[:2] and close(functions._lognorm_w(mu, sigma, x), at[2], 10)


def test_lognorm_ppf_sizes_match_the_tolerance_reason():
    # |y| = |ln x| and sigma |q| = |ln x - mu| at the normal results.
    path = Path(__file__).resolve().parent.parent / "data" / "lognorm_ppf.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()[1:]]
    sizes = [
        (abs(math.log(y)), abs(math.log(y) - r["in"]["mu"][0]), r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    y_max, q_max, p = max(sizes)
    assert p == 5e-324 and 151 < y_max < 152 and 153 < q_max < 154
    assert max(q for _, q, _ in sizes) == q_max
    assert (152 + 154 * 2.4) * 2.0**-53 < 6e-14 < functions.LOGNORM_PPF.tolerance


@pytest.mark.parametrize("alpha", [2.0**-5, 0.5, 1.0, 3.5, 100.0])
@pytest.mark.parametrize("u", [1.0, 1 + 2.0**-52, 1 + 1e-10, 1.5, 1e5, 1e300])
def test_pareto_log_is_the_exponential(alpha, u):
    # log(X / x_m) is exponential with rate alpha: F(x) = expon_cdf(alpha log(x / x_m)).
    # The log is taken at 1200 bits here, not through _mp_pareto_t. The
    # exponential uses the same expm1 / exp, so this checks the log only; the
    # series and quadrature tests below avoid both.
    with mpmath.workdps(50):
        alpha, xm, u = mp(alpha, 2.0**-10, u)
        x = u * xm
        with mpmath.workprec(1200):
            t = alpha * mpmath.log(u)
        assert close(functions._mp_pareto_cdf(alpha, xm, x), functions._mp_expon_cdf(t))
        assert close(functions._mp_pareto_sf(alpha, xm, x), functions._mp_expon_sf(t))


_PARETO_GRID = set(functions.PARETO_CDF.grid())


@pytest.mark.parametrize("alpha", [2.0**-5, 0.1, 1.0, 3.5, 100.0])
@pytest.mark.parametrize("d", [2.0**-52, 1e-15, 1e-10, 1e-5])
def test_pareto_cdf_near_the_scale_matches_the_binomial_series(alpha, d):
    # x = x_m (1 + d): F = 1 - (1 + d)^-alpha = sum_{k>=1} (-1)^(k+1)
    # binom(alpha + k - 1, k) d^k, summed exactly term by term (no log or
    # expm1). Terms fall by a factor of about d (alpha + k) / k, so 40 terms
    # are far below 1e-40 relative for d <= 1e-5 and alpha <= 100.
    with mpmath.workdps(60):
        # The grid's u is the double 1 + d; use its exact distance from 1.
        alpha, d = mpmath.mpf(alpha), mpmath.mpf(1 + d) - 1
        xm = mpmath.mpf(2) ** -10
        x = xm * (1 + d)
        assert (alpha, float(xm), float(x)) in _PARETO_GRID
        series = mpmath.fsum(
            (-1) ** (k + 1) * mpmath.binomial(alpha + k - 1, k) * d**k for k in range(1, 40)
        )
    with mpmath.workdps(50):
        assert close(functions._mp_pareto_cdf(alpha, xm, x), series)


@pytest.mark.parametrize(
    ("alpha", "xm", "x"),
    [(2.0**-5, 2.0**-500, 2.0**-499), (0.1, 1.0, 1 + 1e-10), (0.5, 1024.0, 3000.0),
     (1.5, 2.0**500, 2.0**501), (3.5, 2.0**-10, 0.01), (10.0, 1.0, 1.2), (100.0, 1.0, 1.001)],
)  # fmt: skip
def test_pareto_cdf_and_sf_match_the_integral_of_the_density(alpha, xm, x):
    # F(x) = integral of f(s) = alpha x_m^alpha / s^(alpha + 1) from x_m to x,
    # by quadrature in s (no logarithm or power of the ratio involved).
    with mpmath.workdps(50):
        alpha, xm, x = mp(alpha, xm, x)
        integral = mpmath.quad(
            lambda s: alpha * xm**alpha / s ** (alpha + 1), [xm, (xm + x) / 2, x]
        )
        cdf = functions._mp_pareto_cdf(alpha, xm, x)
        assert close(cdf, integral, digits=35)
        assert close(cdf + functions._mp_pareto_sf(alpha, xm, x), 1)


def test_pareto_far_tail_below_the_smallest_double():
    # x/x_m = 2^1500: Q = 2^(-1500 alpha) exactly, e.g. 2^-46.875 for alpha = 2^-5.
    with mpmath.workdps(50):
        xm, x = mp(2.0**-500, 2.0**1000)
        for alpha in mp(2.0**-5, 0.1, 0.5):
            assert close(functions._mp_pareto_sf(alpha, xm, x), mpmath.mpf(2) ** (-1500 * alpha))


def test_pareto_cdf_and_sf_below_and_at_the_support_start():
    with mpmath.workdps(50):
        xm = mpmath.mpf(2)
        for alpha in mp(0.1, 2.0, 100.0):
            for x in mp(-1e300, 0.0, 1.0, 2 * (1 - 2.0**-53)):
                assert functions._mp_pareto_cdf(alpha, xm, x) == 0
                assert functions._mp_pareto_sf(alpha, xm, x) == 1
            assert functions._mp_pareto_cdf(alpha, xm, xm) == 0
            assert functions._mp_pareto_sf(alpha, xm, xm) == 1


def test_pareto_grid_has_exact_standardised_x():
    # Every x >= x_m is a normal double equal to u x_m with u from _PARETO_U,
    # so x / x_m is exact as a real number (the tolerance reasons rely on it).
    grid = functions.PARETO_CDF.grid()
    assert len(grid) == 1062
    us = {mpmath.ldexp(mantissa, e) for mantissa, e in functions._PARETO_U}
    for alpha, xm, x in grid:
        assert math.frexp(xm)[0] == 0.5  # a power of two
        if x >= xm:
            assert 2.0**-1022 <= x < math.inf
            assert mpmath.mpf(x) / mpmath.mpf(xm) in us
    # The far-tail points are present where x is finite.
    assert (2.0**-5, 2.0**-500, 2.0**1000) in grid
    assert (2.0**-5, 2.0**-10, math.ldexp(1.1, 1020)) in grid


@pytest.mark.parametrize("p", [5e-324, 1e-300, 1e-16, 0.25, 0.5, 0.75, 1 - 1e-10, 1 - 2.0**-53])
@pytest.mark.parametrize(("alpha", "xm"), [(2.0**-5, 2.0**-500), (1.5, 1.0), (100.0, 2.0**500)])
def test_pareto_ppf_inverts_the_cdf_and_matches_the_power(alpha, xm, p):
    # F(ppf(p)) = p (relative to min(p, 1 - p)), and ppf(p) = x_m (1 - p)^(-1/alpha)
    # with 1 - p exact. Both at 1200 bits: for p = 5e-324 the quantile is
    # x_m (1 + p / alpha), so x - x_m needs about 1100 bits to survive.
    with mpmath.workprec(1200):
        alpha, xm, p = mp(alpha, xm, p)
        x = functions._mp_pareto_ppf(alpha, xm, p)
        if p < 0.5:
            assert close(functions._mp_pareto_cdf(alpha, xm, x), p)
        else:
            assert close(functions._mp_pareto_sf(alpha, xm, x), 1 - p)
        route = xm * (1 - p) ** (-1 / alpha)
    with mpmath.workdps(50):
        assert close(functions._mp_pareto_ppf(alpha, xm, p), route)
        assert functions._mp_pareto_ppf(alpha, xm, mpmath.mpf(0)) == xm


def test_pareto_ppf_poles_are_tagged_by_rule():
    assert functions._pareto_ppf_tags(2.0, 1.0, 1.0) == [POLE]
    assert functions._pareto_ppf_tags(2.0, 1.0, 0.0) == []
    assert functions._pareto_ppf_tags(2.0, 1.0, 1 - 2.0**-53) == []


def test_pareto_condition_numbers_match_the_tolerance_reasons():
    # cdf: at most 1 in t; sf: t, below 708.4 at every normal result.
    with mpmath.workdps(30):
        grid = [(a, s, x) for a, s, x in functions.PARETO_CDF.grid() if x > s]
        assert max(functions.pareto_cdf_condition_in_t(*p) for p in grid) <= 1
        normal = [p for p in grid if functions._mp_pareto_sf(*mp(*p)) >= mpmath.mpf(2) ** -1022]
        worst = max(functions.pareto_sf_condition_in_t(*p) for p in normal)
        assert 690 < worst < 708.4
        assert 708.4 * 2.0**-53 < 7.9e-14 < functions.PARETO_SF.tolerance


def _one_over_is_exact(a):
    return Fraction(1) / Fraction(a) == Fraction(1 / a)


def test_pareto_ppf_log_size_matches_the_tolerance_reason():
    # Rounding 1/alpha costs |log(x/x_m)| 2^-53: at most 345.4 2^-53 at a
    # normal result where 1/alpha is not a double (shape 2^-5 has 1/alpha = 32).
    path = Path(__file__).resolve().parent.parent / "data" / "pareto_ppf.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y / r["in"]["xm"][0])), r["in"]["alpha"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags") and not _one_over_is_exact(r["in"]["alpha"][0])
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < math.inf
    ]
    worst, alpha, p = max(logs)
    assert (alpha, p) == (0.1, 1 - 1e-15)
    assert 345.3 < worst < 345.5
    assert worst * 2.0**-53 < 3.9e-14 < functions.PARETO_PPF.tolerance
    assert {a for a in functions._PARETO_SHAPES if _one_over_is_exact(a)} == {
        2.0**-5,
        0.5,
        1.0,
        2.0,
    }


_GUMBEL_GRID = set(functions.GUMBEL_CDF.grid())
_BIG = 2.0**1023


@pytest.mark.parametrize(
    ("mu", "beta", "x"),
    [(0.0, 1.0, -6.5), (0.0, 1.0, -3.0), (0.0, 1.0, 0.0), (0.0, 1.0, 2.5), (0.0, 1.0, 36.75),
     (0.0, 2.0**-500, 2.0**-498), (-10.0, 0.25, -10.5), (100.0, 8.0, 140.0),
     (-_BIG, _BIG, 0.0), (-_BIG, _BIG, _BIG), (-_BIG, _BIG, 1.5 * _BIG)],
)  # fmt: skip
def test_gumbel_cdf_and_sf_match_the_integral_of_the_density(mu, beta, x):
    # F = integral of exp(-t - exp(-t)) dt from -inf to z = (x - mu) / beta, and
    # Q the same from z to inf, by quadrature in t (no closed form involved).
    # The integrand is below exp(-162000) for t < -12 and below exp(-150) times
    # its value at z beyond z + 150, far below 30 digits of F (at least
    # exp(-665) here) and Q.
    with mpmath.workdps(50):
        mu, beta, x = mp(mu, beta, x)
        assert (float(mu), float(beta), float(x)) in _GUMBEL_GRID
        z = (x - mu) / beta

        def density(t):
            return mpmath.exp(-t - mpmath.exp(-t))

        # The integrand is sharp next to z in the lower tail (it falls like
        # exp(-w e^v) at t = z - v, w = exp(-z)), so the breakpoints double
        # away from z on the scale 1 / (1 + w), and away from the peak at 0.
        scale = 1 / (1 + mpmath.exp(-z))
        lower = [z] + [z - 2**k * scale for k in range(80) if z - 2**k * scale > -12] + [-12]
        upper = sorted({z, z + 150} | {t for k in range(-5, 8) for t in (2**k, -(2**k)) if t > z})
        # mpmath.quad stops at an absolute error of about 10^-50, so the
        # integrand is scaled by an exact power of two near 1 / density(z).
        c = mpmath.ldexp(1, -mpmath.mag(density(z)))
        lower = mpmath.quad(lambda t: c * density(t), sorted(lower)) / c
        upper = mpmath.quad(lambda t: c * density(t), upper) / c
        assert close(functions._mp_gumbel_cdf(mu, beta, x), lower, digits=30)
        assert close(functions._mp_gumbel_sf(mu, beta, x), upper, digits=30)


@pytest.mark.parametrize("z", [5.0, 36.75, 700.0, 745.0, 1e5])
def test_gumbel_sf_upper_tail_matches_the_series(z):
    # Q = 1 - exp(-w) = sum_{k>=1} (-1)^(k+1) w^k / k! with w = exp(-z) <= 0.0068,
    # summed term by term (no expm1); 40 terms are far below 1e-40 relative.
    with mpmath.workdps(60):
        w = mpmath.exp(-mpmath.mpf(z))
        series = mpmath.fsum((-1) ** (k + 1) * w**k / mpmath.factorial(k) for k in range(1, 40))
    with mpmath.workdps(50):
        assert (0.0, 1.0, z) in _GUMBEL_GRID
        assert close(functions._mp_gumbel_sf(*mp(0.0, 1.0, z)), series)


def test_gumbel_is_the_exponential_of_exp_minus_z():
    # exp(-Z) is exponential with rate 1, so F(x) = expon_sf(exp(-z)) and
    # Q(x) = expon_cdf(exp(-z)). The exponential uses the same exp / expm1, so
    # this checks the standardisation only.
    with mpmath.workdps(50):
        for mu, beta, x in [(0.0, 1.0, -2.0), (-10.0, 0.25, -9.0), (-_BIG, _BIG, _BIG)]:
            mu, beta, x = mp(mu, beta, x)
            with mpmath.workprec(1200):
                w = mpmath.exp(-(x - mu) / beta)
            assert close(functions._mp_gumbel_cdf(mu, beta, x), functions._mp_expon_sf(w))
            assert close(functions._mp_gumbel_sf(mu, beta, x), functions._mp_expon_cdf(w))


def test_gumbel_grid_has_exact_standardised_x():
    # Every x is mu + beta z with beta a power of two and z from _GUMBEL_Z, so
    # z = (x - mu) / beta is exact (the tolerance reasons rely on it).
    grid = functions.GUMBEL_CDF.grid()
    assert len(grid) == 1092
    zs = {mpmath.mpf(z) for z in functions._GUMBEL_Z}
    with mpmath.workprec(2200):
        for mu, beta, x in grid:
            assert math.frexp(beta)[0] == 0.5
            assert x == 0 or 2.0**-1022 <= abs(x) < math.inf
            assert (mpmath.mpf(x) - mu) / beta in zs
    # The points where x - mu exceeds the largest double although z is 2 to 2.75.
    big = [x for mu, _, x in grid if mu == -_BIG and x - _BIG / 2 >= _BIG / 2]
    assert big == [_BIG, 1.25 * _BIG, 1.5 * _BIG, 1.75 * _BIG]


def test_gumbel_ppf_matches_the_root_of_the_cdf():
    # The quantile solves F(x) = p: found here by the secant method at 1200 bits
    # from the cdf (exp only, no logarithm), starting next to the route's value.
    for mu, beta in [(0.0, 1.0), (-10.0, 0.25), (-_BIG, _BIG)]:
        for p in (5e-324, 1e-300, 0.25, 0.9, 0.95, 1 - 1e-15, *functions._GUMBEL_NEAR_INV_E):
            with mpmath.workdps(50):
                route = functions._mp_gumbel_ppf(*mp(mu, beta, p))
            with mpmath.workprec(1200):
                m, b, q = mp(mu, beta, p)
                start = route * (1 + mpmath.mpf(10) ** -30) + b * mpmath.mpf(10) ** -30
                root = mpmath.findroot(
                    lambda x, m=m, b=b, q=q: functions._mp_gumbel_cdf(m, b, x) - q, (route, start)
                )
            with mpmath.workdps(50):
                assert close(route, root)


def test_gumbel_ppf_poles_are_tagged_by_rule():
    assert functions._gumbel_ppf_tags(0.0, 1.0, 0.0) == [POLE]
    assert functions._gumbel_ppf_tags(0.0, 1.0, 1.0) == [POLE]
    assert functions._gumbel_ppf_tags(0.0, 1.0, 5e-324) == []
    assert functions._gumbel_ppf_tags(0.0, 1.0, 1 - 2.0**-53) == []


def test_gumbel_condition_numbers_match_the_tolerance_reasons():
    # cdf: w = exp(-z) in w, below 708.4 at any normal result (665.1 at z = -6.5
    # on this grid); sf: at most 1.
    with mpmath.workdps(30):
        grid = functions.GUMBEL_CDF.grid()
        normal = [p for p in grid if functions._mp_gumbel_cdf(*mp(*p)) >= mpmath.mpf(2) ** -1022]
        worst = max(functions.gumbel_cdf_condition_in_w(*p) for p in normal)
        assert 665 < worst < 665.2
        assert 708.4 * 2.0**-53 < 7.9e-14 < functions.GUMBEL_CDF.tolerance
        assert max(functions.gumbel_sf_condition_in_w(*p) for p in grid) <= 1


def test_gumbel_ppf_cancellation_matches_the_tolerance_reason():
    # mu - beta L costs (|mu| + |beta L|) / |mu - beta L| ulps: at most about 14
    # at a finite, untagged quantile on this grid.
    with mpmath.workdps(30):
        ratios = []
        for mu, beta, p in functions.GUMBEL_PPF.grid():
            if functions._gumbel_ppf_tags(mu, beta, p):
                continue
            m, b, q = mp(mu, beta, p)
            beta_l = b * mpmath.log(-mpmath.log(q))
            if abs(m - beta_l) < mpmath.mpf(2) ** 1024:
                ratios.append(((abs(m) + abs(beta_l)) / abs(m - beta_l), mu, beta, p))
        worst, mu, beta, p = max(ratios)
        assert (mu, beta, p) == (-10.0, 0.25, 1 - 1e-15)
        assert 13 < worst < 14
        assert 14 * 2.0**-53 < 1.6e-15 < functions.GUMBEL_PPF.tolerance


def test_gumbel_ppf_near_inverse_e_cancels_within_the_extra_precision():
    # log(-log p) loses the leading bits of 1 + log p near p = 1/e: |1 + log p|
    # is at least 1e-17 on this grid, so about 57 of the 60 extra bits.
    with mpmath.workprec(300):
        worst = min(abs(1 + mpmath.log(mpmath.mpf(p))) for p in functions._GUMBEL_NEAR_INV_E)
    assert 2.0**-60 < 1e-17 < worst


_LAPLACE_GRID = set(functions.LAPLACE_CDF.grid())


@pytest.mark.parametrize(
    ("mu", "b", "x"),
    [(0.0, 1.0, -740.0), (0.0, 1.0, -3.0), (0.0, 1.0, -(2.0**-60)), (0.0, 1.0, 0.0),
     (0.0, 1.0, 2.0**-30), (0.0, 1.0, 2.5), (0.0, 1.0, 36.75), (0.0, 1.0, 700.0),
     (0.0, 2.0**-500, 2.0**-498), (-10.0, 0.25, -10.5), (100.0, 8.0, 140.0),
     (-_BIG, _BIG, -0.75 * _BIG), (-_BIG, _BIG, _BIG), (-_BIG, _BIG, 1.75 * _BIG)],
)  # fmt: skip
def test_laplace_cdf_and_sf_match_the_integral_of_the_density(mu, b, x):
    # F = integral of e^-|t| / 2 dt from -inf to z = (x - mu) / b and Q the
    # same from z to inf, by quadrature in t (no closed form involved), split
    # at 0, where the density has a kink. Below z - 200 (or above z + 200) the
    # integrand is below e^-200 times its value at z.
    with mpmath.workdps(50):
        mu, b, x = mp(mu, b, x)
        assert (float(mu), float(b), float(x)) in _LAPLACE_GRID
        z = (x - mu) / b

        def density(t):
            return mpmath.exp(-abs(t)) / 2

        # mpmath.quad stops at an absolute error of about 10^-50, so the
        # integrand is scaled by an exact power of two near 1 / density(z).
        c = mpmath.ldexp(1, -mpmath.mag(density(z)))

        def integral(lo, hi):
            points = sorted({lo, hi} | ({mpmath.mpf(0)} if lo < 0 < hi else set()))
            return mpmath.quad(lambda t: c * density(t), points) / c

        lower = integral(min(z, 0) - 200, z) if z <= 0 else 1 - integral(z, max(z, 0) + 200)
        upper = integral(z, max(z, 0) + 200) if z >= 0 else 1 - integral(min(z, 0) - 200, z)
        assert close(functions._mp_laplace_cdf(mu, b, x), lower, digits=30)
        assert close(functions._mp_laplace_sf(mu, b, x), upper, digits=30)


def test_laplace_cdf_tail_is_half_the_exponential_sf():
    # For z >= 0, Q = expon_sf(z) / 2 and F(-z) = Q(z) (symmetry); checks the
    # branches and the mirror in _mp_laplace_sf, not exp itself.
    with mpmath.workdps(50):
        for z in (0.0, 2.0**-60, 3.0, 700.0, 1e5):
            q = functions._mp_expon_sf(mpmath.mpf(z)) / 2
            assert close(functions._mp_laplace_sf(*mp(0.0, 1.0, z)), q)
            assert close(functions._mp_laplace_cdf(*mp(0.0, 1.0, -z)), q)
            assert close(functions._mp_laplace_cdf(*mp(0.0, 1.0, z)), 1 - q)


def test_laplace_grid_has_exact_standardised_x():
    # Every x is mu + b z with b a power of two and z from _LAPLACE_Z, so
    # z = (x - mu) / b is exact (the tolerance reasons rely on it).
    grid = functions.LAPLACE_CDF.grid()
    assert len(grid) == 1828
    zs = {mpmath.mpf(z) for z in functions._LAPLACE_Z}
    with mpmath.workprec(2200):
        for mu, b, x in grid:
            assert math.frexp(b)[0] == 0.5
            assert x == 0 or 2.0**-1022 <= abs(x) < math.inf
            assert (mpmath.mpf(x) - mu) / b in zs
    # The points where x - mu exceeds the largest double although z is 2 to 2.75.
    big = [x for mu, _, x in grid if mu == -_BIG and x - _BIG / 2 >= _BIG / 2]
    assert big == [_BIG, 1.25 * _BIG, 1.5 * _BIG, 1.75 * _BIG]


def test_laplace_ppf_matches_the_root_of_the_cdf():
    # The quantile solves F(x) = p: found here by the secant method at 1200 bits
    # from the cdf (exp only, no logarithm), starting next to the route's value.
    for mu, b in [(0.0, 1.0), (100.0, 8.0), (-_BIG, _BIG)]:
        for p in (5e-324, 1e-300, 1e-5, 0.25, 0.5 - 2.0**-54, 0.5, 0.5 + 2.0**-53, 0.95, 1 - 1e-15):
            with mpmath.workdps(50):
                route = functions._mp_laplace_ppf(*mp(mu, b, p))
            with mpmath.workprec(1200):
                m, bb, q = mp(mu, b, p)
                start = route * (1 + mpmath.mpf(10) ** -30) + bb * mpmath.mpf(10) ** -30
                root = mpmath.findroot(
                    lambda x, m=m, bb=bb, q=q: functions._mp_laplace_cdf(m, bb, x) - q,
                    (route, start),
                )
            with mpmath.workdps(50):
                assert close(route, root)


def test_laplace_ppf_poles_are_tagged_by_rule():
    assert functions._laplace_ppf_tags(0.0, 1.0, 0.0) == [POLE]
    assert functions._laplace_ppf_tags(0.0, 1.0, 1.0) == [POLE]
    assert functions._laplace_ppf_tags(0.0, 1.0, 5e-324) == []
    assert functions._laplace_ppf_tags(0.0, 1.0, 1 - 2.0**-53) == []


def test_laplace_ppf_cancellation_matches_the_tolerance_reason():
    # mu + b L costs (|mu| + |b L|) / |mu + b L| ulps: at most about 14 at a
    # finite, nonzero, untagged quantile on this grid.
    with mpmath.workdps(30):
        ratios = []
        for mu, b, p in functions.LAPLACE_PPF.grid():
            if functions._laplace_ppf_tags(mu, b, p):
                continue
            x = functions._mp_laplace_ppf(*mp(mu, b, p))
            if x != 0 and abs(x) < mpmath.mpf(2) ** 1024:
                ratios.append((functions.laplace_ppf_cancellation(mu, b, p), mu, b, p))
        worst, mu, b, p = max(ratios)
        assert (mu, b, p) == (100.0, 8.0, 1e-5)
        assert 13 < worst < 14
        assert 14 * 2.0**-53 < 1.6e-15 < functions.LAPLACE_PPF.tolerance


_CAUCHY_GRID = set(functions.CAUCHY_CDF.grid())


@pytest.mark.parametrize(
    ("mu", "b", "x"),
    [(0.0, 1.0, -1e308), (0.0, 1.0, -1e16), (0.0, 1.0, -40.0), (0.0, 1.0, -3.0),
     (0.0, 1.0, -(2.0**-60)), (0.0, 1.0, 0.0), (0.0, 1.0, 2.0**-30), (0.0, 1.0, 2.5),
     (0.0, 1.0, 1e5), (0.0, 2.0**-500, -(2.0**-500) * 1e305), (-10.0, 0.25, -10.5),
     (100.0, 8.0, 140.0), (-_BIG, _BIG, -0.75 * _BIG), (-_BIG, _BIG, _BIG),
     (-_BIG, _BIG, 1.75 * _BIG)],
)  # fmt: skip
def test_cauchy_cdf_and_sf_match_the_integral_of_the_density(mu, b, x):
    # F = integral of 1 / (pi (1 + t^2)) dt from -inf to z = (x - mu) / b and Q
    # the same from z to inf, by quadrature (no arctangent involved). The
    # substitution t = 1 / u maps the tail beyond |z| >= 1 to (0, 1 / |z|], where
    # the integrand 1 / (pi (1 + u^2)) is smooth; the rest is integrated in t.
    with mpmath.workdps(50):
        mu, b, x = mp(mu, b, x)
        assert (float(mu), float(b), float(x)) in _CAUCHY_GRID
        z = (x - mu) / b

        def density(t):
            return 1 / (mpmath.pi * (1 + t**2))

        def tail(a):
            # integral from a >= 0 to inf.
            if a >= 1:
                # u = s / a, so the integral over s in [0, 1] is near 1 / pi
                # (quad stops at an absolute error of about 10^-50).
                return mpmath.quad(lambda s: density(s / a), [0, 1]) / a
            return mpmath.quad(density, [a, 1]) + mpmath.quad(density, [0, 1])

        lower = tail(-z) if z <= 0 else 1 - tail(z)
        upper = tail(z) if z >= 0 else 1 - tail(-z)
        assert close(functions._mp_cauchy_cdf(mu, b, x), lower, digits=30)
        assert close(functions._mp_cauchy_sf(mu, b, x), upper, digits=30)


def test_cauchy_cdf_tail_matches_its_series():
    # For z > 1, Q(z) = F(-z) = (1 / pi) (1/z - 1/(3 z^3) + 1/(5 z^5) - ...),
    # the series of acot; checks the branches and the mirror, not acot itself.
    with mpmath.workdps(50):
        for z in (1e3, 1e16, 1e100, 1e308):
            zz = mpmath.mpf(z)
            q = sum((-1) ** k / ((2 * k + 1) * zz ** (2 * k + 1)) for k in range(40))
            q /= mpmath.pi
            assert close(functions._mp_cauchy_sf(*mp(0.0, 1.0, z)), q)
            assert close(functions._mp_cauchy_cdf(*mp(0.0, 1.0, -z)), q)
            assert close(functions._mp_cauchy_cdf(*mp(0.0, 1.0, z)), 1 - q)


def test_cauchy_grid_has_exact_standardised_x():
    # Every x is mu + b z with b a power of two and z from _CAUCHY_Z, so
    # z = (x - mu) / b is exact (the tolerance reasons rely on it).
    grid = functions.CAUCHY_CDF.grid()
    assert len(grid) == 1762
    zs = {mpmath.mpf(z) for z in functions._CAUCHY_Z}
    with mpmath.workprec(2200):
        for mu, b, x in grid:
            assert math.frexp(b)[0] == 0.5
            assert x == 0 or 2.0**-1022 <= abs(x) < math.inf
            assert (mpmath.mpf(x) - mu) / b in zs
    # The points where x - mu exceeds the largest double although z is 2 to 2.75.
    big = [x for mu, _, x in grid if mu == -_BIG and x - _BIG / 2 >= _BIG / 2]
    assert big == [_BIG, 1.25 * _BIG, 1.5 * _BIG, 1.75 * _BIG]


def test_cauchy_ppf_matches_the_root_of_the_cdf():
    # The quantile solves F(x) = p: found here by the secant method at 1200 bits
    # from the cdf (acot, no tangent), starting next to the route's value.
    for mu, b in [(0.0, 1.0), (100.0, 8.0), (-_BIG, _BIG)]:
        for p in (5e-324, 1e-300, 1e-5, 0.125, 0.5 - 2.0**-54, 0.5 + 2.0**-53, 0.875, 1 - 1e-15):
            with mpmath.workdps(50):
                route = functions._mp_cauchy_ppf(*mp(mu, b, p))
            with mpmath.workprec(1200):
                m, bb, q = mp(mu, b, p)
                start = route * (1 + mpmath.mpf(10) ** -30) + bb * mpmath.mpf(10) ** -30
                root = mpmath.findroot(
                    lambda x, m=m, bb=bb, q=q: functions._mp_cauchy_cdf(m, bb, x) - q,
                    (route, start),
                )
            with mpmath.workdps(50):
                assert close(route, root)
    # Exactly the location at p = 1/2.
    assert functions._mp_cauchy_ppf(*mp(100.0, 8.0, 0.5)) == 100


def test_cauchy_ppf_poles_are_tagged_by_rule():
    assert functions._cauchy_ppf_tags(0.0, 1.0, 0.0) == [POLE]
    assert functions._cauchy_ppf_tags(0.0, 1.0, 1.0) == [POLE]
    assert functions._cauchy_ppf_tags(0.0, 1.0, 5e-324) == []
    assert functions._cauchy_ppf_tags(0.0, 1.0, 1 - 2.0**-53) == []


def test_cauchy_ppf_cancellation_matches_the_tolerance_reason():
    # mu + b T costs (|mu| + |b T|) / |mu + b T| ulps: at most about 8.8 at a
    # finite, nonzero, untagged quantile on this grid. The quantile is 0 only
    # where mu = 0, p = 1/2 (T = 0) and mu = -b, p = 3/4 (T = 1).
    with mpmath.workdps(30):
        ratios, zeros = [], []
        for mu, b, p in functions.CAUCHY_PPF.grid():
            if functions._cauchy_ppf_tags(mu, b, p):
                continue
            x = functions._mp_cauchy_ppf(*mp(mu, b, p))
            if x == 0:
                zeros.append((mu, b, p))
            elif abs(x) < mpmath.mpf(2) ** 1024:
                ratios.append((functions.cauchy_ppf_cancellation(mu, b, p), mu, b, p))
        assert sorted(zeros) == [
            (-_BIG, _BIG, 0.75), (0.0, 2.0**-500, 0.5), (0.0, 1.0, 0.5), (0.0, 2.0**500, 0.5)
        ]  # fmt: skip
        worst, mu, b, p = max(ratios)
        assert (mu, b, p) == (-10.0, 0.25, 0.99)
        assert 8.7 < worst < 8.8
        assert 8.8 * 2.0**-53 < 1e-15 < functions.CAUCHY_PPF.tolerance


_LEVY_GRID = set(functions.LEVY_CDF.grid())


@pytest.mark.parametrize(
    ("mu", "c", "x"),
    [(0.0, 1.0, 46 / 65536), (0.0, 1.0, 5 / 4096), (0.0, 1.0, 0.25), (0.0, 1.0, 1.0),
     (0.0, 1.0, 40.0), (0.0, 1.0, 1e10), (0.0, 1.0, 1e308), (0.0, 2.0**-500, 2.0**1023),
     (0.0, 2.0**500, 2.0**500 / 64), (-10.0, 0.25, -9.5), (100.0, 8.0, 140.0),
     (100.0, 8.0, 100 + 8 * 56 / 65536), (-_BIG, _BIG, -_BIG + _BIG / 1024),
     (-_BIG, _BIG, _BIG / 4), (-_BIG, _BIG, _BIG), (-_BIG, _BIG, 1.75 * _BIG)],
)  # fmt: skip
def test_levy_cdf_and_sf_match_the_integral_of_the_density(mu, c, x):
    # The density in t = (x - mu) / c is (2 pi)^-1/2 t^-3/2 e^(-1 / (2 tau)) at tau.
    # Substituting tau = t / r^2 turns Q = integral from t to inf into
    # 2 (2 pi t)^-1/2 times the integral of e^(-r^2 / (2 t)) over r in [0, 1], and F =
    # integral from 0 to t into the same over r in [1, inf); with r = 1 + v the latter
    # is e^(-1 / (2 t)) times the integral of e^(-(2 v + v^2) / (2 t)) over v >= 0, so
    # the tiny lower tail keeps its relative precision. No error function involved.
    with mpmath.workdps(50):
        mu, c, x = mp(mu, c, x)
        assert (float(mu), float(c), float(x)) in _LEVY_GRID
        t = (x - mu) / c
        k = 2 / mpmath.sqrt(2 * mpmath.pi * t)
        w = min(t, 1)  # the width of the peaks of both integrands
        upper = k * mpmath.quad(lambda r: mpmath.exp(-(r**2) / (2 * t)), [0, w / 2, 1])
        if t >= 1:
            # F >= erfc(1 / sqrt 2) = 0.317, so 1 - Q does not cancel (and the v
            # integrand would spread out to v of about sqrt(t)).
            lower = 1 - upper
        else:
            v_integral = mpmath.quad(
                lambda v: mpmath.exp(-(2 * v + v**2) / (2 * t)), [0, w, 1, mpmath.inf]
            )
            lower = k * mpmath.exp(-1 / (2 * t)) * v_integral
        assert close(functions._mp_levy_cdf(mu, c, x), lower, digits=30)
        assert close(functions._mp_levy_sf(mu, c, x), upper, digits=30)


def test_levy_cdf_and_sf_outside_the_support():
    # F = 0 and Q = 1 for x <= mu: grid points with t = -1, -2^-30 and 0.
    for mu, c, x in [
        (0.0, 1.0, -1.0),
        (0.0, 1.0, -(2.0**-30)),
        (0.0, 1.0, 0.0),
        (100.0, 8.0, 100.0),
    ]:
        assert (mu, c, x) in _LEVY_GRID
        assert functions._mp_levy_cdf(*mp(mu, c, x)) == 0
        assert functions._mp_levy_sf(*mp(mu, c, x)) == 1


def test_levy_sf_upper_tail_matches_its_series():
    # For large t, Q = erf(s) with s = (2 t)^-1/2 small: Q = (2 / sqrt(pi)) sum
    # (-1)^k s^(2k+1) / (k! (2k+1)); checks the far tail where c / (2 (x - mu)) is
    # subnormal (3 2^550, 3 2^570) or below every double, and F = 1 - Q there.
    with mpmath.workdps(50):
        for mu, c, x in [
            (0.0, 1.0, 1e308),
            (0.0, 2.0**-500, 3 * 2.0**550),
            (0.0, 2.0**-500, 3 * 2.0**570),
            (0.0, 2.0**-500, 2.0**700),
            (0.0, 2.0**-500, 2.0**1023),
        ]:
            assert (mu, c, x) in _LEVY_GRID
            s = 1 / mpmath.sqrt(2 * (mpmath.mpf(x) - mu) / c)
            q = sum(
                (-1) ** k * s ** (2 * k + 1) / (mpmath.factorial(k) * (2 * k + 1))
                for k in range(30)
            )
            q *= 2 / mpmath.sqrt(mpmath.pi)
            assert close(functions._mp_levy_sf(*mp(mu, c, x)), q)
            assert close(functions._mp_levy_cdf(*mp(mu, c, x)), 1 - q)


def test_levy_grid_has_exact_standardised_x():
    # Every x is mu + c t with c a power of two and t from _LEVY_T, so
    # t = (x - mu) / c is exact (the tolerance reasons rely on it).
    grid = functions.LEVY_CDF.grid()
    assert len(grid) == 1776
    ts = set(functions._LEVY_T)
    with mpmath.workprec(2200):
        for mu, c, x in grid:
            assert math.frexp(c)[0] == 0.5
            assert x == 0 or 2.0**-1022 <= abs(x) < math.inf
            assert (mpmath.mpf(x) - mu) / c in ts
    # The points where x - mu exceeds the largest double although t is 2 to 2.75.
    big = [x for mu, _, x in grid if mu == -_BIG and x >= _BIG]
    assert big == [_BIG, 1.25 * _BIG, 1.5 * _BIG, 1.75 * _BIG]
    # The scale 2^-500 reaches t = 3 2^1030 to 2^1523, where c / (2 (x - mu)) is below
    # the smallest normal double (and from 2^1100 below every double).
    far = [x for mu, c, x in grid if c == 2.0**-500 and x >= 2.0**530]
    assert far == [3 * 2.0**530, 3 * 2.0**550, 3 * 2.0**570] + [
        2.0**k for k in (600, 700, 1000, 1023)
    ]


def test_levy_cdf_condition_matches_the_tolerance_reason():
    # The largest condition number in u at a normal F below 1 is about 698 (as the
    # tolerance reason says), and Q's is at most 1/2.
    with mpmath.workdps(30):
        conds = []
        for mu, c, x in functions.LEVY_CDF.grid():
            if x <= mu:
                continue
            f = functions._mp_levy_cdf(*mp(mu, c, x))
            if mpmath.mpf(2) ** -1022 <= f < 1 - mpmath.mpf(2) ** -53:
                conds.append((functions.levy_cdf_condition_in_u(mu, c, x), mu, c, x))
            # At most 1/2 (its limit for tiny u), up to the rounding of this sum.
            assert functions.levy_sf_condition_in_u(mu, c, x) <= 0.5 + mpmath.mpf(10) ** -25
        worst, mu, c, x = max(conds)
        assert (mu, c, x) == (100.0, 8.0, 100 + 8 * 47 / 65536)
        assert 697 < worst < 698
        assert 698 * 2.0**-53 < 7.8e-14 < functions.LEVY_CDF.tolerance


def test_levy_ppf_matches_the_root_of_the_cdf():
    # The quantile solves F(x) = p: found here by the secant method at 1200 bits
    # from the cdf (erfc, no inverse), starting next to the route's value.
    for mu, c in [(0.0, 1.0), (100.0, 8.0), (-_BIG, _BIG)]:
        for p in (5e-324, 1e-300, 1e-5, 0.25, 0.5 - 2.0**-54, 0.5 + 2.0**-53, 0.9, 1 - 2.0**-53):
            with mpmath.workdps(50):
                route = functions._mp_levy_ppf(*mp(mu, c, p))
            with mpmath.workprec(1200):
                m, cc, q = mp(mu, c, p)
                start = route * (1 + mpmath.mpf(10) ** -30) + cc * mpmath.mpf(10) ** -30
                root = mpmath.findroot(
                    lambda x, m=m, cc=cc, q=q: functions._mp_levy_cdf(m, cc, x) - q,
                    (route, start),
                )
            with mpmath.workdps(50):
                assert close(route, root)
    # Exactly the location at p = 0.
    assert functions._mp_levy_ppf(*mp(100.0, 8.0, 0.0)) == 100


def test_levy_ppf_poles_are_tagged_by_rule():
    assert functions._levy_ppf_tags(0.0, 1.0, 1.0) == [POLE]
    assert functions._levy_ppf_tags(0.0, 1.0, 0.0) == []
    assert functions._levy_ppf_tags(0.0, 1.0, 5e-324) == []
    assert functions._levy_ppf_tags(0.0, 1.0, 1 - 2.0**-53) == []


def test_levy_ppf_cancellation_matches_the_tolerance_reason():
    # mu + c / z^2 costs (|mu| + |c / z^2|) / |mu + c / z^2| ulps: at most about 7.2 at
    # a finite, nonzero, untagged quantile on this grid. It is 0 only at p = 0 with
    # mu = 0.
    with mpmath.workdps(30):
        ratios, zeros = [], []
        for mu, c, p in functions.LEVY_PPF.grid():
            if functions._levy_ppf_tags(mu, c, p):
                continue
            x = functions._mp_levy_ppf(*mp(mu, c, p))
            if x == 0:
                zeros.append((mu, c, p))
            elif abs(x) < mpmath.mpf(2) ** 1024:
                ratios.append((functions.levy_ppf_cancellation(mu, c, p), mu, c, p))
        assert sorted(zeros) == [(0.0, 2.0**-500, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0**500, 0.0)]
        worst, mu, c, p = max(ratios)
        assert (mu, c, p) == (-_BIG, _BIG, 0.25)
        assert 7.1 < worst < 7.2
        assert 7.2 * 2.0**-53 < 8e-16 < functions.LEVY_PPF.tolerance


@pytest.mark.parametrize(
    ("a", "b", "p"),
    [
        (2.0, 3.0, 0.25),
        (0.1, 0.1, 1e-20),
        (0.5, 10.0, 0.999),
        (2.0, 10.0, 1e-300),  # SciPy returns NaN
        (1000.0, 0.5, 5e-324),  # SciPy returns 0.375
        (100.0, 10.0, 1e-300),  # SciPy is off by 1.6e-10
        (1000.0, 10000.0, 0.999),  # SciPy is off by 2.8e-13
        (10000.0, 10000.0, 0.25),
        (10000.0, 0.1, 1 - 1e-15),
        (0.1, 1000.0, 1 - 1e-10),
    ],
)
def test_betaincinv_inverts_the_series(a, b, p):
    # I_x(a, b) = p at the quantile, with I from the positive-term series
    # (_mp_betainc), not the continued fraction the quantile search uses; in
    # the upper half 1 - I = 1 - p, with 1 - I from the mirrored series so
    # that it keeps its relative precision. 250 digits so that 1 - x is
    # accurate even where it is about 1e-150 (b = 0.1 in the upper tail).
    with mpmath.workdps(250):
        a, b, p = mp(a, b, p)
        x = functions._mp_betaincinv(a, b, p)
        assert 0 < x < 1
        if p <= 0.5:
            assert close(functions._mp_betainc(a, b, x), p, digits=35)
        else:
            assert close(functions._mp_betainc(b, a, 1 - x), 1 - p, digits=35)


@pytest.mark.parametrize("p", [1e-300, 1e-10, 0.25, 0.5, 0.75, 1 - 1e-15])
def test_betaincinv_matches_closed_forms(p):
    # b = 1: I = x^a, x = p^(1/a); a = 1: I = 1 - (1 - x)^b, x = 1 - (1 - p)^(1/b);
    # a = b = 1/2: I = (2/pi) asin(sqrt x), x = sin(pi p / 2)^2.
    with mpmath.workdps(60):
        (p,) = mp(p)
        three, half = mpmath.mpf(3), mpmath.mpf(1) / 2
        inv = functions._mp_betaincinv
        assert close(inv(three, mpmath.mpf(1), p), p ** (1 / three), digits=50)
        assert close(inv(mpmath.mpf(1), three, p), -mpmath.expm1(mpmath.log1p(-p) / 3), digits=50)
        assert close(inv(half, half, p), mpmath.sin(mpmath.pi * p / 2) ** 2, digits=50)


def test_betaincinv_ends_of_the_support_are_not_poles():
    # p = 0 and p = 1 give the ends of the support, 0 and 1: finite, so untagged.
    assert functions._mp_betaincinv(mpmath.mpf(2), mpmath.mpf(3), mpmath.mpf(0)) == 0
    assert functions._mp_betaincinv(mpmath.mpf(2), mpmath.mpf(3), mpmath.mpf(1)) == 1
    assert functions.BETAINCINV.input_tags(2.0, 3.0, 0.0) == []
    assert functions.BETAINCINV.input_tags(2.0, 3.0, 1.0) == []


def test_betaincinv_log_magnitude_matches_the_tolerance_reason():
    # The largest |log x| at a nonzero, normal quantile below 1 in the committed
    # data (a = 1, b = 1e4, p = 1e-300, x = 1e-304), against the 2 |L| 2^-53 bound.
    data = Path(__file__).resolve().parent.parent / "data" / "betaincinv.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    logs = [
        (abs(math.log(y)), r["in"]["a"][0], r["in"]["b"][0], r["in"]["p"][0])
        for r in rows
        if not r.get("tags")
        for y in [float(r["mpmath"]["y"][0])]
        if 2.0**-1022 <= y < 1
    ]
    worst, a, b, p = max(logs)
    assert (a, b, p) == (1.0, 1e4, 1e-300)
    assert 699.9 < worst < 700.0
    assert functions.BETAINCINV.tolerance >= 2 * worst * 2.0**-53


# --- beta_cdf, beta_sf ---


@pytest.mark.parametrize(
    ("a", "b", "x"),
    [(0.5, 2.0, 1e-5), (2.0, 5.0, 0.3), (5.0, 2.0, 0.95), (10.0, 100.0, 0.2), (0.1, 0.5, 0.999)],
)
def test_beta_cdf_and_sf_match_mpmath_betainc(a, b, x):
    # Both sides of the switch point against mpmath's hypergeometric
    # evaluation, which converges quickly for these moderate parameters.
    with mpmath.workdps(50):
        a, b, x = mp(a, b, x)
        assert close(functions._mp_beta_cdf(a, b, x), mpmath.betainc(a, b, 0, x, regularized=True))
        assert close(functions._mp_beta_sf(a, b, x), mpmath.betainc(a, b, x, 1, regularized=True))


@pytest.mark.parametrize(
    ("a", "b", "x"),
    [
        (1e4, 1e4, 0.46464554479100484),
        (1e4, 0.1, 0.9998951379361433),
        (0.1, 1e4, 0.00010486206385669913),
        (1e3, 1e4, 0.06350021450445349),
        (100.0, 0.1, 1 - 1e-15),
        (10.0, 1e4, 1e-20),
    ],
)
def test_beta_cdf_and_sf_match_the_series_in_both_tails(a, b, x):
    # The positive-term series DLMF 8.17.8 (_mp_betainc) on the smaller tail,
    # at the points where statrs fails worst: I directly, 1 - I as the
    # mirrored series I_(1-x)(b, a), so neither forms 1 minus a number near 1.
    with mpmath.workdps(50):
        a, b, x = mp(a, b, x)
        y = mpmath.fsub(1, x, exact=True)
        assert close(functions._mp_beta_cdf(a, b, x), functions._mp_betainc(a, b, x))
        assert close(functions._mp_beta_sf(a, b, x), functions._mp_betainc(b, a, y))


@pytest.mark.parametrize("x", [1e-300, 1e-20, 1e-5, 0.3, 0.5, 0.9, 1 - 1e-10, 1 - 1e-15])
def test_beta_cdf_and_sf_match_closed_forms(x):
    # a = b = 1/2: I = (2/pi) asin(sqrt x), 1 - I = (2/pi) asin(sqrt(1 - x))
    # (where SciPy is wrong: x = 1e-20 for sf, 1 - 1e-15 for cdf);
    # b = 1: I = x^a; a = 1: 1 - I = (1 - x)^b.
    with mpmath.workdps(60):
        (x,) = mp(x)
        y = mpmath.fsub(1, x, exact=True)
        half, three = mpmath.mpf(1) / 2, mpmath.mpf(3)
        one = mpmath.mpf(1)
        assert close(
            functions._mp_beta_cdf(half, half, x), 2 / mpmath.pi * mpmath.asin(mpmath.sqrt(x))
        )
        assert close(
            functions._mp_beta_sf(half, half, x), 2 / mpmath.pi * mpmath.asin(mpmath.sqrt(y))
        )
        assert close(functions._mp_beta_cdf(three, one, x), x**three)
        assert close(functions._mp_beta_sf(one, three, x), y**three)


def test_beta_cdf_and_sf_sum_to_one_and_cover_the_ends():
    with mpmath.workdps(50):
        for a, b in (mp(0.1, 1e4), mp(2.0, 3.0), mp(1e4, 1e4)):
            assert functions._mp_beta_cdf(a, b, mpmath.mpf(0)) == 0
            assert functions._mp_beta_cdf(a, b, mpmath.mpf(1)) == 1
            assert functions._mp_beta_sf(a, b, mpmath.mpf(0)) == 1
            assert functions._mp_beta_sf(a, b, mpmath.mpf(1)) == 0
            for x in mp(1e-300, 0.01, 0.5, 0.99, 1 - 1e-15):
                total = functions._mp_beta_cdf(a, b, x) + functions._mp_beta_sf(a, b, x)
                assert close(total, 1)


def test_beta_grid_covers_both_tails_the_ends_and_the_bulk():
    grid = functions._beta_x_grid()
    assert all(0 <= x <= 1 for _, _, x in grid)
    assert len(grid) == len(set(grid)) == 1895
    narrow = [x for a, b, x in grid if (a, b) == (1e4, 1e4)]
    # The mean +- k sd points of the narrowest distribution are on the grid.
    assert sum(0.46 < x < 0.54 for x in narrow) == 7
    assert {0.0, 1.0, 1e-300, 1 - 1e-15} <= set(narrow)


def test_beta_prefix_logarithm_matches_the_tolerance_reason():
    # Over the committed data: the prefix belongs to the smaller of I and
    # 1 - I, so its rounding error, relative to the dataset's value v, is the
    # largest prefix logarithm times min(I, 1 - I) / v (times 2^-53). The
    # worst is ln Gamma(2e4) = 1.78e5 at a = b = 1e4, in both datasets.
    data = Path(__file__).resolve().parent.parent / "data"
    rows = {}
    for name in ("beta_cdf", "beta_sf"):
        lines = (data / f"{name}.jsonl").read_text().splitlines()[1:]
        rows[name] = [json.loads(line) for line in lines]
    for name in ("beta_cdf", "beta_sf"):
        worst = (0.0, None)
        for rc, rs in zip(rows["beta_cdf"], rows["beta_sf"], strict=True):
            a, b, x = (rc["in"][k][0] for k in "abx")
            i, q = float(rc["mpmath"]["y"][0]), float(rs["mpmath"]["y"][0])
            v = i if name == "beta_cdf" else q
            if 0 < x < 1 and v >= 2.0**-1022:
                loss = functions.beta_prefix_log_magnitude(a, b, x) * min(i, q) / v
                worst = max(worst, (loss, (a, b)))
        loss, at = worst
        assert at == (1e4, 1e4)
        assert 1.78e5 < loss < 1.781e5
        assert 1.97e-11 < loss * 2.0**-53 < 1.98e-11  # "2.0e-11"
        assert loss * 2.0**-53 > functions.BETA_CDF.tolerance == functions.BETA_SF.tolerance


# --- f_cdf, f_sf, f_ppf ---


@pytest.mark.parametrize(
    ("d1", "d2", "x"), [(1.0, 2.0, 0.5), (5.0, 10.0, 2.0), (10.0, 5.0, 0.3), (30.0, 100.0, 1.5)]
)
def test_f_cdf_and_sf_match_mpmath_betainc(d1, d2, x):
    # Both sides of the switch point against mpmath's hypergeometric
    # evaluation of I_z(d1/2, d2/2) with z = d1 x / (d1 x + d2).
    with mpmath.workdps(50):
        d1, d2, x = mp(d1, d2, x)
        z = d1 * x / (d1 * x + d2)
        a, b = d1 / 2, d2 / 2
        assert close(functions._mp_f_cdf(d1, d2, x), mpmath.betainc(a, b, 0, z, regularized=True))
        assert close(functions._mp_f_sf(d1, d2, x), mpmath.betainc(a, b, z, 1, regularized=True))


@pytest.mark.parametrize(
    ("d1", "d2", "x"),
    [
        (1e5, 10.0, 0.5),
        (1e5, 30.0, 0.1),
        (30.0, 1e5, 3.5823761667632135),
        (1e5, 1e5, 1.0),
        (1e5, 0.5, 1e10),
        (0.5, 1e5, 1e-10),
    ],
)
def test_f_cdf_and_sf_match_the_series_on_the_smaller_tail(d1, d2, x):
    # The positive-term series DLMF 8.17.8 (_mp_betainc) on the smaller tail,
    # at points where SciPy fails (rows 1 to 3) or statrs fails worst: I_z(a, b)
    # for the lower tail, I_y(b, a) for the upper tail, y = d2 / (d1 x + d2).
    with mpmath.workdps(50):
        d1, d2, x = mp(d1, d2, x)
        a, b = d1 / 2, d2 / 2
        z, y = d1 * x / (d1 * x + d2), d2 / (d1 * x + d2)
        cdf, sf = functions._mp_f_cdf(d1, d2, x), functions._mp_f_sf(d1, d2, x)
        if cdf <= sf:
            assert close(cdf, functions._mp_betainc(a, b, z))
        else:
            assert close(sf, functions._mp_betainc(b, a, y))


@pytest.mark.parametrize("x", [1e-300, 1e-20, 1e-5, 0.5, 1.0, 3.0, 1e5, 1e20, 1e300])
def test_f_cdf_and_sf_match_closed_forms(x):
    # d1 = d2 = 1: F = (2/pi) atan(sqrt x), 1 - F = (2/pi) atan(1 / sqrt x)
    # (SciPy is wrong at x = 1e20 for cdf and 1e-20 for sf);
    # d2 = 2: F = z^(d1/2); d1 = 2: 1 - F = (1 + 2 x / d2)^(-d2/2).
    with mpmath.workdps(60):
        (x,) = mp(x)
        one, two, five = mp(1, 2, 5)
        root = mpmath.sqrt(x)
        assert close(functions._mp_f_cdf(one, one, x), 2 / mpmath.pi * mpmath.atan(root))
        assert close(functions._mp_f_sf(one, one, x), 2 / mpmath.pi * mpmath.atan(1 / root))
        z = five * x / (five * x + two)
        assert close(functions._mp_f_cdf(five, two, x), z ** (five / 2))
        assert close(functions._mp_f_sf(two, five, x), (1 + 2 * x / five) ** (-five / 2))


def test_f_cdf_and_sf_sum_to_one_and_cover_zero():
    with mpmath.workdps(50):
        for d1, d2 in (mp(0.5, 1e5), mp(2.0, 3.0), mp(1e5, 1e5), mp(1e5, 0.5)):
            assert functions._mp_f_cdf(d1, d2, mpmath.mpf(0)) == 0
            assert functions._mp_f_sf(d1, d2, mpmath.mpf(0)) == 1
            for x in mp(1e-300, 0.01, 1.0, 50.0, 1e300):
                total = functions._mp_f_cdf(d1, d2, x) + functions._mp_f_sf(d1, d2, x)
                assert close(total, 1)


def test_f_grid_covers_both_tails_and_the_bulk():
    grid = functions._f_x_grid()
    assert all(x >= 0 for _, _, x in grid)
    assert len(grid) == len(set(grid)) == 1999
    narrow = [x for d1, d2, x in grid if (d1, d2) == (1e5, 1e5)]
    # The 1 + k s points of the narrowest distribution (s = 0.0063) are on the grid.
    assert sum(0.9 < x < 1.1 and x != 1 for x in narrow) == 6
    assert {0.0, 1e-300, 1.0, 1e300} <= set(narrow)


def _f_rows(name):
    data = Path(__file__).resolve().parent.parent / "data"
    lines = (data / f"{name}.jsonl").read_text().splitlines()[1:]
    return [json.loads(line) for line in lines]


@pytest.mark.parametrize(("name", "upper", "worst", "at"), [
    ("f_cdf", False, 1647.79, (1e5, 1e5, 0.9367544467966324)),
    ("f_sf", True, 1548.63, (1e5, 1e5, 1.0632455532033676)),
])  # fmt: skip
def test_f_condition_matches_the_tolerance_reason(name, upper, worst, at):
    # Over the committed data: the largest |x f / F| (or over 1 - F) at a
    # normal result, which a rounding of z or y multiplies.
    best = (0.0, None)
    for row in _f_rows(name):
        d1, d2, x = (row["in"][k][0] for k in ("d1", "d2", "x"))
        if x > 0 and float(row["mpmath"]["y"][0]) >= 2.0**-1022:
            best = max(best, (float(functions.f_condition(d1, d2, x, upper)), (d1, d2, x)))
    assert best[1] == at
    assert worst < best[0] < worst + 0.01
    assert 2 * best[0] * 2.0**-53 < functions.F_CDF.tolerance == functions.F_SF.tolerance


def test_f_prefix_logarithm_matches_the_tolerance_reason():
    # The largest prefix logarithm: ln Gamma(a + b) at d1 = d2 = 1e5, 1.05e6,
    # whose rounding alone (1.2e-10) exceeds the tolerance.
    worst = functions.beta_prefix_log_magnitude(5e4, 5e4, 0.5)
    assert worst == pytest.approx(float(mpmath.loggamma(1e5)))
    assert 1.05e6 < worst < 1.06e6
    assert 1.1e-10 < worst * 2.0**-53 < 1.2e-10
    assert worst * 2.0**-53 > functions.F_CDF.tolerance


@pytest.mark.parametrize(
    "p", [5e-324, 1e-300, 1e-20, 1e-5, 0.1, 0.5, 0.9, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15]
)
def test_f_ppf_matches_closed_forms(p):
    # d1 = d2 = 1: x = tan(pi p / 2)^2; d1 = 2: x = (d2 / 2) ((1 - p)^(-2/d2) - 1);
    # d2 = 2: z = p^(2/d1), x = (2 / d1) z / (1 - z).
    with mpmath.workdps(60):
        (p,) = mp(p)
        one, two, five = mp(1, 2, 5)
        assert close(functions._mp_f_ppf(one, one, p), mpmath.tan(mpmath.pi * p / 2) ** 2)
        q = mpmath.fsub(1, p, exact=True)
        assert close(
            functions._mp_f_ppf(two, five, p), five / 2 * mpmath.expm1(-2 / five * mpmath.log(q))
        )
        z = p ** (2 / five)
        assert close(
            functions._mp_f_ppf(five, two, p),
            2 / five * z / -mpmath.expm1(2 / five * mpmath.log(p)),
        )


@pytest.mark.parametrize(
    ("d1", "d2", "p"),
    [
        (1e5, 1.0, 5e-324),
        (1000.0, 1e5, 5e-324),
        (10.0, 30.0, 5e-324),
        (1e5, 0.5, 0.5),
        (0.5, 1e5, 1 - 1e-15),
        (1e5, 1e5, 0.999),
        (5.0, 0.5, 1e-300),
        (1000.0, 0.5, 1e-100),
        (0.5, 1e5, 1e-100),
    ],
)
def test_f_ppf_inverts_the_cdf(d1, d2, p):
    # F(ppf(p)) = p, through the independently tested cdf and sf, at points
    # where SciPy fails (rows 1 to 3 and the last, whose quantile is far below
    # the smallest double) or statrs fails worst.
    with mpmath.workdps(50):
        d1, d2, p = mp(d1, d2, p)
        x = functions._mp_f_ppf(d1, d2, p)
        if p <= 0.5:
            assert close(functions._mp_f_cdf(d1, d2, x), p, 35)
        else:
            assert close(functions._mp_f_sf(d1, d2, x), 1 - p, 35)


def test_f_ppf_ends_and_tags():
    assert functions._f_ppf_tags(1.0, 1.0, 1.0) == [POLE]
    assert functions._f_ppf_tags(1.0, 1.0, 0.0) == []
    with mpmath.workdps(50):
        assert functions._mp_f_ppf(*mp(2.0, 3.0, 0.0)) == 0
    assert len(functions._f_ppf_grid()) == 81 * 20


def test_f_ppf_log_magnitude_matches_the_tolerance_reason():
    worst = (0.0, None)
    for row in _f_rows("f_ppf"):
        v = row["mpmath"]["y"][0]
        if not isinstance(v, str) and 2.0**-1022 <= v:
            worst = max(
                worst, (abs(math.log(v)), tuple(row["in"][k][0] for k in ("d1", "d2", "p")))
            )
    assert worst[1][0] == 2.0 and worst[1][2] == 1e-300
    assert 690.7 < worst[0] < 690.8
    assert functions.F_PPF.tolerance >= 2 * worst[0] * 2.0**-53


# --- gammaln, gammasgn, gammainc_large ---


@pytest.mark.parametrize("x", [0.25, 0.75, 3.5, 30.25, 170.5, 1e-300])
def test_gammaln_matches_log_of_gamma_for_positive_x(x):
    with mpmath.workdps(50):
        (x,) = mp(x)
        assert close(functions._mp_gammaln(x), mpmath.log(mpmath.gamma(x)))


@pytest.mark.parametrize("x", [-1e-300, -0.5, -2.75, -2.5, -29.75, -1000000.5])
def test_gammaln_matches_reflection_formula_for_negative_x(x):
    # log|Gamma(x)| = log(pi / |sin(pi x)|) - log Gamma(1 - x) (DLMF 5.5.3).
    with mpmath.workdps(60):
        (x,) = mp(x)
        expected = mpmath.log(mpmath.pi / abs(mpmath.sinpi(x))) - mpmath.loggamma(1 - x)
        assert close(functions._mp_gammaln(x), expected, digits=45)


@pytest.mark.parametrize("d", [1e-5, -1e-5, 1e-10, -1e-10, 1e-15, -1e-15])
def test_gammaln_matches_the_series_next_to_one_and_two(d):
    # log Gamma(1 + z) = -gamma z + sum_{k>=2} (-1)^k zeta(k) z^k / k (DLMF 5.7.3),
    # and log Gamma(2 + z) = log Gamma(1 + z) + log(1 + z). z = x - 1 or x - 2 is
    # exact for the double x. SciPy is up to 29% off at these points.
    with mpmath.workdps(50):
        for zero in (1.0, 2.0):
            x = mpmath.mpf(zero + d)
            z = x - zero
            series = -mpmath.euler * z + mpmath.nsum(
                lambda k, z=z: (-1) ** k * mpmath.zeta(k) * z**k / k, [2, mpmath.inf]
            )
            if zero == 2.0:
                series += mpmath.log1p(z)
            assert close(functions._mp_gammaln(x), series)


@pytest.mark.parametrize("x", [1e10, 1e100, 1e300, 2.5e305, 1.7976931348623157e308])
def test_gammaln_matches_stirling_at_large_x(x):
    # (x - 1/2) log x - x + log(2 pi) / 2 + 1 / (12 x) - 1 / (360 x^3) (DLMF 5.11.1);
    # the next term is below 1e-40 for x >= 1e10.
    with mpmath.workdps(60):
        (x,) = mp(x)
        stirling = (x - 0.5) * mpmath.log(x) - x + mpmath.log(2 * mpmath.pi) / 2
        stirling += 1 / (12 * x) - 1 / (360 * x**3)
        assert close(functions._mp_gammaln(x), stirling, digits=45)


def test_gammaln_reflection_amplification_matches_the_tolerance_reason():
    # The largest cancellation in the reflection formula on the grid (x = -2.75,
    # next to the zero at -2.7477), against the tolerance's 2.5x margin.
    with mpmath.workdps(30):
        worst, at = max(
            (functions.gammaln_reflection_amplification(x), x)
            for (x,) in functions.GAMMALN.grid()
            if x < 0 and x != int(x)
        )
        zero = mpmath.findroot(lambda t: functions._mp_gammaln(t), -2.75)
    assert at == -2.75
    assert 663 < worst < 664
    assert -2.7477 < zero < -2.7476
    assert 7.3e-14 < worst * 2.0**-53 < 7.4e-14
    assert functions.GAMMALN.tolerance >= 2.5 * worst * 2.0**-53


def test_gammaln_overflows_where_the_grid_description_says():
    # log Gamma(x) = 1.797e308 at x of about 2.553e305.
    data = Path(__file__).resolve().parent.parent / "data" / "gammaln.jsonl"
    rows = [json.loads(line) for line in data.read_text().splitlines()[1:]]
    overflow = sorted(r["in"]["x"][0] for r in rows if "overflow" in r.get("tags", []))
    assert overflow == [2.6e305, 1e306, 1e307, 1e308, 1.7976931348623157e308]
    with mpmath.workdps(30):
        assert functions._mp_gammaln(mpmath.mpf(2.5e305)) < mpmath.mpf(1.7976931348623157e308)
        assert functions._mp_gammaln(mpmath.mpf(2.6e305)) > mpmath.mpf("1.8e308")


def test_gammaln_poles_are_tagged_by_rule_and_left_out_of_gammasgn():
    assert functions.GAMMALN.input_tags(0.0) == [POLE]
    assert functions.GAMMALN.input_tags(-30.0) == [POLE]
    assert functions.GAMMALN.input_tags(-0.25) == []
    assert functions.GAMMALN.input_tags(1.0) == []
    poles = {(x,) for (x,) in functions.GAMMALN.grid() if functions.GAMMALN.input_tags(x)}
    assert len(poles) == 31
    assert set(functions.GAMMASGN.grid()) == set(functions.GAMMALN.grid()) - poles


def test_gammasgn_matches_the_parity_rule():
    # +1 for x > 0; (-1)^ceil(-x) for negative non-integers (DLMF 5.2.1, 5.5.3).
    with mpmath.workdps(50):
        for (x,) in functions.GAMMASGN.grid():
            expected = 1 if x > 0 else (-1) ** math.ceil(-x)
            assert functions._mp_gammasgn(mpmath.mpf(x)) == expected


def _quad_gammainc(a: float, x: float, lower: bool) -> mpmath.mpf:
    # P(a, x) = x^a e^-x / Gamma(a) int_0^1 exp((a - 1) log(1 - s) + x s) ds (x < a - 1), or
    # Q(a, x) = x^a e^-x / Gamma(a) int_0^inf exp((a - 1) log(1 + s) - x s) ds, with
    # t = x (1 -+ s) in DLMF 8.2.1 / 8.2.2. The integrand peaks at s = 0 and decays
    # on the scale h, so the breakpoints double from h / 8.
    a, x = mp(a, x)
    with mpmath.extraprec(80):
        h = 1 / (abs(a - 1 - x) + mpmath.sqrt(a))
        if lower:
            ends = [h * 2**j / 8 for j in range(200) if h * 2**j / 8 < 1]
            points = [0, *ends, 1]
            integral = mpmath.quad(lambda s: mpmath.exp((a - 1) * mpmath.log1p(-s) + x * s), points)
        else:
            points = [0, *[h * 2**j / 8 for j in range(40)], mpmath.inf]
            integral = mpmath.quad(lambda s: mpmath.exp((a - 1) * mpmath.log1p(s) - x * s), points)
        return mpmath.exp(a * mpmath.log(x) - x - mpmath.loggamma(a)) * integral


@pytest.mark.parametrize(
    ("a", "k"),
    [(1e8, -30), (1e8, -10), (1e8, -3), (1e9, -30), (1e9, -10), (1e9, -1), (1e8, 0), (1e9, 1)],
)
def test_gammainc_large_matches_numerical_integration(a, k):
    # An independent route at the bulk points, including the four where SciPy
    # is off (k = -10 and -30). Q for x >= a - 1, where P = 1 - Q is not small.
    x = a + k * math.sqrt(a)
    with mpmath.workdps(30):
        a_mp, x_mp = mp(a, x)
        if x < a - 1:
            assert close(functions._mp_gamma_cdf(a_mp, x_mp), _quad_gammainc(a, x, True), 25)
        else:
            assert close(functions._mp_gamma_sf(a_mp, x_mp), _quad_gammainc(a, x, False), 25)


def test_gammainc_large_condition_number_matches_the_tolerance_reason():
    # The largest condition number at a normal, nonzero point is 9.5e5 at a = 1e9,
    # x = a - 30 sqrt(a): rounding x / a there would cost about 1e-10, far above
    # the tolerance, while a careful method's k^2 2^-53 fits.
    with mpmath.workdps(30):
        worst, at = max(
            (functions.gamma_cdf_condition(a, x), (a, x))
            for a, x in functions.GAMMAINC_LARGE.grid()
            if functions._mp_gamma_cdf(mpmath.mpf(a), mpmath.mpf(x)) >= mpmath.mpf(2) ** -1022
        )
    assert at == (1e9, 1e9 - 30 * math.sqrt(1e9))
    assert 9.4e5 < worst < 9.6e5
    assert 1.0e-10 < worst * 2.0**-53 < 1.1e-10
    assert 30**2 * 2.0**-53 < functions.GAMMAINC_LARGE.tolerance / 9


def test_gammainc_large_grid_crosses_1e8_and_reaches_1e20():
    grid = functions.GAMMAINC_LARGE.grid()
    assert min(a for a, _ in grid) == 1e8
    assert max(a for a, _ in grid) == 1e20
    # Points on both sides of x = 1e8 at a = 1e8, and x > 1e8 below the mean.
    assert any(a == 1e8 and x < 1e8 for a, x in grid)
    assert any(a == 1e8 and x > 1e8 for a, x in grid)
    assert all(any(a == b and 1e8 < x < a for b, x in grid) for a in (1e9, 1e10, 1e15, 1e20))


# --- digamma, trigamma ---


def _harmonic(n: int) -> mpmath.mpf:
    return mpmath.fsum(mpmath.mpf(1) / k for k in range(1, n + 1))


@pytest.mark.parametrize("n", [1, 2, 7, 30, 180])
def test_digamma_and_trigamma_match_closed_forms_at_integers_and_half_integers(n):
    # psi(n) = H_{n-1} - gamma, psi(n + 1/2) = -gamma - 2 log 2 + sum_{k=1}^{n} 2 / (2k - 1),
    # psi'(n) = pi^2 / 6 - sum_{k<n} 1 / k^2, psi'(n + 1/2) = pi^2 / 2 - sum 4 / (2k - 1)^2
    # (DLMF 5.4.14, 5.4.15, 5.15.2 with 5.15.5).
    with mpmath.workdps(50):
        half = mpmath.mpf(n) + mpmath.mpf(0.5)
        odd = [mpmath.mpf(2 * k - 1) for k in range(1, n + 1)]
        assert close(mpmath.digamma(n), _harmonic(n - 1) - mpmath.euler)
        expected = -mpmath.euler - 2 * mpmath.log(2) + mpmath.fsum(2 / o for o in odd)
        assert close(mpmath.digamma(half), expected)
        expected = mpmath.pi**2 / 6 - mpmath.fsum(mpmath.mpf(1) / k**2 for k in range(1, n))
        assert close(functions._mp_trigamma(mpmath.mpf(n)), expected)
        expected = mpmath.pi**2 / 2 - mpmath.fsum(4 / o**2 for o in odd)
        assert close(functions._mp_trigamma(half), expected)


def _digamma_by_recurrence(x: mpmath.mpf, n: int) -> mpmath.mpf:
    # psi(x) = psi(x + n) - sum_{k=0}^{n-1} 1 / (x + k) (DLMF 5.5.2), evaluated by
    # mpmath at x + n > 0, away from where the dataset's point lies.
    return mpmath.digamma(x + n) - mpmath.fsum(1 / (x + k) for k in range(n))


def _trigamma_by_recurrence(x: mpmath.mpf, n: int) -> mpmath.mpf:
    return functions._mp_trigamma(x + n) + mpmath.fsum(1 / (x + k) ** 2 for k in range(n))


@pytest.mark.parametrize("x", [-0.25, -2.5, -22.75, -29.75, -100.5, -1000.25, -1e-300])
def test_digamma_and_trigamma_match_the_recurrence_for_negative_x(x):
    with mpmath.workdps(60):
        (x,) = mp(x)
        n = int(-mpmath.floor(x)) + 20
        assert close(mpmath.digamma(x), _digamma_by_recurrence(x, n), digits=45)
        assert close(functions._mp_trigamma(x), _trigamma_by_recurrence(x, n), digits=45)


@pytest.mark.parametrize("x", [-10000.5, -1000000.5, -999999999999999.5])
def test_digamma_matches_the_reflection_formula_at_large_negative_x(x):
    # psi(x) = psi(1 - x) - pi cot(pi x) (DLMF 5.5.4). The recurrence would need
    # up to 1e15 terms here.
    with mpmath.workdps(60):
        (x,) = mp(x)
        cot = mpmath.cospi(x) / mpmath.sinpi(x)
        assert close(mpmath.digamma(x), mpmath.digamma(1 - x) - mpmath.pi * cot, digits=45)


@pytest.mark.parametrize("x", [-0.25, -2.5, -29.75, -100.5, -10000.5])
def test_trigamma_reflection_matches_mpmaths_own_route_for_negative_x(x):
    # The dataset's reflection route against mpmath.psi(1, x), which recurs up
    # from x and is too slow to use below about -1e4.
    with mpmath.workdps(50):
        (x,) = mp(x)
        assert close(functions._mp_trigamma(x), mpmath.psi(1, x))


def test_digamma_matches_the_recurrence_next_to_its_zero():
    # Next to the zero the value is tiny (2e-17 at the double nearest x0); the
    # recurrence from x + 40 subtracts terms near 3.7 and keeps 40 of 60 digits.
    with mpmath.workdps(60):
        for (x,) in functions.DIGAMMA.grid():
            if abs(x - functions._DIGAMMA_ROOT) < 1e-4:
                x = mpmath.mpf(x)
                assert close(mpmath.digamma(x), _digamma_by_recurrence(x, 40), digits=38)
        root = mpmath.findroot(mpmath.digamma, 1.46)
    assert float(root) == functions._DIGAMMA_ROOT


@pytest.mark.parametrize("x", [1e10, 1e100, 1e300, 1.7976931348623157e308])
def test_digamma_and_trigamma_match_the_asymptotic_series_at_large_x(x):
    # psi(x) = log x - 1 / (2x) - 1 / (12 x^2) + 1 / (120 x^4) - ...,
    # psi'(x) = 1 / x + 1 / (2 x^2) + 1 / (6 x^3) - 1 / (30 x^5) + ... (DLMF 5.11.2,
    # 5.15.8); the next terms are below 1e-40 relative for x >= 1e10.
    with mpmath.workdps(60):
        (x,) = mp(x)
        psi = mpmath.log(x) - 1 / (2 * x) - 1 / (12 * x**2) + 1 / (120 * x**4)
        assert close(mpmath.digamma(x), psi, digits=45)
        psi1 = 1 / x + 1 / (2 * x**2) + 1 / (6 * x**3) - 1 / (30 * x**5)
        assert close(functions._mp_trigamma(x), psi1, digits=45)


def test_digamma_reflection_amplification_matches_the_tolerance_reason():
    # The largest cancellation in the reflection formula on the grid (x = -22.75,
    # next to the zero at -22.7502), against the tolerance's 2.5x margin.
    with mpmath.workdps(30):
        worst, at = max(
            (functions.digamma_reflection_amplification(x), x)
            for (x,) in functions.DIGAMMA.grid()
            if x < 0 and x != int(x)
        )
        zero = mpmath.findroot(mpmath.digamma, -22.75)
    assert at == -22.75
    assert 1312 < worst < 1314
    assert -22.7503 < zero < -22.7502
    assert 1.4e-13 < worst * 2.0**-53 < 1.5e-13
    assert functions.DIGAMMA.tolerance >= 2.5 * worst * 2.0**-53


def test_trigamma_reflection_amplification_is_small():
    # (pi^2 / sin^2 + psi'(1 - x)) / psi'(x) <= (pi^2 + pi^2 / 6) / (pi^2 - pi^2 / 6) = 1.4,
    # since psi'(1 - x) <= psi'(1) = pi^2 / 6 for x <= 0.
    with mpmath.workdps(30):
        for (x,) in functions.TRIGAMMA.grid():
            if x < 0 and x != int(x) and x > -1e-150:
                x = mpmath.mpf(x)
                t1 = mpmath.pi**2 / mpmath.sinpi(x) ** 2
                t2 = functions._mp_trigamma(1 - x)
                assert (t1 + t2) / (t1 - t2) <= 1.4


def test_polygamma_grids_and_poles():
    digamma = set(functions.DIGAMMA.grid())
    trigamma = set(functions.TRIGAMMA.grid())
    near_root = digamma - trigamma - {(-999999999999999.5,)}
    assert len(near_root) == 5
    assert all(abs(x - 1.4616321449683623) <= 1.000001e-5 for (x,) in near_root)
    poles = {(x,) for (x,) in trigamma if functions.TRIGAMMA.input_tags(x)}
    assert poles == {(-float(n),) for n in range(31)}
    assert functions.DIGAMMA.input_tags(-3.0) == [POLE]
    assert (-999999999999999.5,) in digamma
    assert min(trigamma) == (-1000000.5,)


def _exact_t_squared(x, mu):
    # (t^2, mean - mu) with t^2 = n (mean - mu)^2 / s^2 in exact rational
    # arithmetic (the doubles are exact binary fractions), a route that shares
    # no code with mpmath's.
    from fractions import Fraction

    xs = [Fraction(v) for v in x]
    n = len(xs)
    mean = sum(xs) / n
    var = sum((v - mean) ** 2 for v in xs) / (n - 1)
    diff = mean - Fraction(mu)
    return n * diff**2 / var, diff


@pytest.mark.parametrize("point", functions._ttest_grid()[::7])
def test_ttest_statistic_matches_exact_rational_arithmetic(point):
    x, mu, alt = point
    t2, diff = _exact_t_squared(x, mu)
    with mpmath.workdps(50):
        t, _ = functions._mp_ttest_1samp(mp(*x), mpmath.mpf(mu), mpmath.mpf(alt))
        exact = mpmath.mpf(t2.numerator) / t2.denominator
        assert close(t * t, exact)
        assert (t > 0) == (diff > 0)


@pytest.mark.parametrize("point", [p for p in functions._ttest_grid() if len(p[0]) in (2, 3)])
def test_ttest_pvalue_matches_closed_forms_for_df_one_and_two(point):
    # df = 1: sf(t) = 1/2 - atan(t) / pi; df = 2: sf(t) = 1/2 - t / (2 sqrt(t^2 + 2)),
    # written as 1 / (sqrt(t^2 + 2) (sqrt(t^2 + 2) + t)) to keep tiny tails exact.
    x, mu, alt = point
    with mpmath.workdps(50):
        t, p = functions._mp_ttest_1samp(mp(*x), mpmath.mpf(mu), mpmath.mpf(alt))

        def sf(u):
            if len(x) == 2:
                return (
                    mpmath.atan(1 / u) / mpmath.pi if u > 0 else 1 - mpmath.atan(-1 / u) / mpmath.pi
                )
            r = mpmath.sqrt(u * u + 2)
            return 1 / (r * (r + u)) if u > 0 else 1 - 1 / (r * (r - u))

        expected = {0.0: 2 * sf(abs(t)), 1.0: sf(t), -1.0: sf(-t)}[alt]
        assert close(p, expected, digits=35)


def test_ttest_rounding_bound_matches_the_tolerance_reason():
    with mpmath.workdps(30):
        worst, point = max((functions.ttest_rounding_bound(*p), p) for p in functions._ttest_grid())
        kappa_t, _ = functions.ttest_condition(*point)
        assert 4.4e-12 < worst <= functions.TTEST_1SAMP.tolerance
        assert (len(point[0]), point[1]) == (100, 9.99)
        assert 3.9e4 < kappa_t < 4.1e4
        assert max(functions.ttest_condition(*p)[1] for p in functions._ttest_grid()) < 99.1


def _exact_chisquare(obs, exp):
    from fractions import Fraction

    return sum(
        (Fraction(o) - Fraction(e)) ** 2 / Fraction(e) for o, e in zip(obs, exp, strict=True)
    )


@pytest.mark.parametrize("point", functions._chisquare_grid()[::5])
def test_chisquare_statistic_matches_exact_rational_arithmetic(point):
    obs, exp, ddof = point
    exact = _exact_chisquare(obs, exp)
    with mpmath.workdps(50):
        stat, _ = functions._mp_chisquare(mp(*obs), mp(*exp), mpmath.mpf(ddof))
        assert close(stat, mpmath.mpf(exact.numerator) / exact.denominator)


@pytest.mark.parametrize(
    "point", [p for p in functions._chisquare_grid() if len(p[0]) - 1 - p[2] in (1, 2)]
)
def test_chisquare_pvalue_matches_closed_forms_and_gammainc(point):
    # dof = 1: Q = erfc(sqrt(x / 2)); dof = 2: Q = e^(-x / 2); any dof:
    # mpmath.gammainc(dof / 2, x / 2, inf).
    obs, exp, ddof = point
    dof = len(obs) - 1 - ddof
    with mpmath.workdps(50):
        stat, p = functions._mp_chisquare(mp(*obs), mp(*exp), mpmath.mpf(ddof))
        closed = mpmath.erfc(mpmath.sqrt(stat / 2)) if dof == 1 else mpmath.exp(-stat / 2)
        assert close(p, closed)
        assert close(
            p, mpmath.gammainc(mpmath.mpf(dof) / 2, stat / 2, mpmath.inf, regularized=True)
        )


def test_chisquare_rounding_bound_matches_the_tolerance_reason():
    with mpmath.workdps(30):
        bounds = []
        for point in functions._chisquare_grid():
            _, p = functions._mp_chisquare(mp(*point[0]), mp(*point[1]), mpmath.mpf(point[2]))
            if p >= mpmath.mpf(2) ** -1022:  # p-values that round to 0 are not at risk
                bounds.append((functions.chisquare_rounding_bound(*point), point))
        worst, point = max(bounds)
    assert 4.2e-12 < worst <= functions.CHISQUARE.tolerance
    assert len(point[0]) == 30


def _exact_anova_f(x, sizes):
    from fractions import Fraction

    groups = functions._anova_groups([Fraction(v) for v in x], sizes)
    n, k = len(x), len(groups)
    mean = sum(sum(g) for g in groups) / n
    means = [sum(g) / len(g) for g in groups]
    ssb = sum(len(g) * (m - mean) ** 2 for g, m in zip(groups, means, strict=True))
    ssw = sum((v - m) ** 2 for g, m in zip(groups, means, strict=True) for v in g)
    return (ssb / (k - 1)) / (ssw / (n - k))


@pytest.mark.parametrize("point", functions._f_oneway_grid()[::3])
def test_f_oneway_statistic_matches_exact_rational_arithmetic(point):
    exact = _exact_anova_f(*point)
    with mpmath.workdps(50):
        f, _ = functions._mp_f_oneway(mp(*point[0]), mp(*point[1]))
        assert close(f, mpmath.mpf(exact.numerator) / exact.denominator)


@pytest.mark.parametrize("point", functions._f_oneway_grid())
def test_f_oneway_pvalue_matches_betainc_and_the_two_sample_t_test(point):
    # p = I_y(d2 / 2, d1 / 2) with y = d2 / (d2 + d1 F) by mpmath.betainc, which
    # includes the points where SciPy is off (p below 1e-52); for two groups
    # F = t^2 for the pooled two-sample t, so p = 2 sf_t(sqrt F) with N - 2
    # degrees of freedom.
    x, sizes = point
    with mpmath.workdps(50):
        f, p = functions._mp_f_oneway(mp(*x), mp(*sizes))
        d1, d2 = mpmath.mpf(len(sizes) - 1), mpmath.mpf(len(x) - len(sizes))
        with mpmath.extraprec(60):
            y = d2 / (d2 + d1 * f)
            beta = mpmath.betainc(d2 / 2, d1 / 2, 0, y, regularized=True)
        assert close(p, beta, digits=35)
        if len(sizes) == 2:
            assert close(p, 2 * functions._mp_t_sf(d2, mpmath.sqrt(f)))


def test_f_oneway_rounding_bound_matches_the_tolerance_reason():
    with mpmath.workdps(30):
        worst, point = max(
            (functions.f_oneway_rounding_bound(*p), p) for p in functions._f_oneway_grid()
        )
    assert 4.1e-12 < worst <= functions.F_ONEWAY.tolerance
    assert point[1] == (30.0, 30.0, 30.0) and point[0][0] == 99.5


@pytest.mark.parametrize("x", [s for s in functions._skewtest_samples()[::3] if len(s) <= 50])
def test_skewness_matches_exact_rational_arithmetic(x):
    # b1^2 = m3^2 / m2^3 with the biased central moments, exactly.
    from fractions import Fraction

    xs = [Fraction(v) for v in x]
    n = len(xs)
    mean = sum(xs) / n
    m2 = sum((v - mean) ** 2 for v in xs) / n
    m3 = sum((v - mean) ** 3 for v in xs) / n
    with mpmath.workdps(50):
        b1 = functions._mp_skew_root_b1(mp(*x))
        exact = m3**2 / m2**3
        assert close(b1 * b1, mpmath.mpf(exact.numerator) / exact.denominator) or exact == 0
        assert (b1 > 0) == (m3 > 0) and (b1 == 0) == (m3 == 0)


def test_skewtest_matches_the_published_example():
    # SciPy's documented example (also in statrs's tests): z = 2.7788579769903414,
    # two-sided p = 0.005455036974740185, both to all digits shown.
    x = [148.0, 154.0, 158.0, 160.0, 161.0, 162.0, 166.0, 170.0, 182.0, 195.0, 236.0]
    with mpmath.workdps(50):
        z, p = functions._mp_skewtest(mp(*x), mpmath.mpf(0))
        assert abs(z - mpmath.mpf("2.7788579769903414")) < 1e-15
        assert abs(p - mpmath.mpf("0.005455036974740185")) < 1e-17


@pytest.mark.parametrize("alt", [-1.0, 0.0, 1.0])
def test_skewtest_pvalue_matches_the_normal_cdf(alt):
    # The erfc forms against mpmath.ncdf (a different function) on both sides.
    with mpmath.workdps(50):
        for z in mp(-12.1, -3.4, -0.03, 0.09, 2.6, 12.1):
            phi = mpmath.ncdf
            expected = {-1.0: phi(z), 1.0: phi(-z), 0.0: 2 * phi(-abs(z))}[alt]
            assert close(functions._mp_norm_tail(z, alt), expected, digits=35)


def test_skewtest_symmetric_samples_have_zero_statistic_and_exact_double_moments():
    # The symmetric samples (values loc + i - (n + 1) / 2) have z = 0, p = 1 (two-sided)
    # and 1/2; and every operation in their moments is exact in double arithmetic
    # (half-integers, cubes far below 2^53), so a careful method returns z = 0 too.
    symmetric = [x for x in functions._skewtest_samples() if x[1] - x[0] == 1.0]
    assert len(symmetric) == 2 * len(functions._SKEW_N)
    for x in symmetric:
        with mpmath.workdps(50):
            assert functions._mp_skewtest(mp(*x), mpmath.mpf(0)) == (0, 1)
        mean = sum(x) / len(x)
        assert mean * len(x) == math.fsum(x)  # the double mean is exact
        assert sum((v - mean) ** 3 for v in x) == 0.0


def test_skewtest_rounding_bound_matches_the_tolerance_reason():
    with mpmath.workdps(30):
        worst, point = max(
            (functions.skewtest_rounding_bound(*p), p) for p in functions._skewtest_grid()
        )
    assert 1.6e-12 < worst <= functions.SKEWTEST.tolerance
    assert len(point[0]) == 100


@pytest.mark.parametrize("n1, n2", [(1, 1), (3, 6), (4, 5), (5, 4), (6, 3), (8, 8), (2, 11)])
def test_mannwhitneyu_null_counts_match_enumeration(n1, n2):
    # The recurrence against every split of the ranks 1..n1+n2 (U = R1 - n1 (n1 + 1) / 2).
    counts = [0] * (n1 * n2 + 1)
    for split in itertools.combinations(range(1, n1 + n2 + 1), n1):
        counts[sum(split) - n1 * (n1 + 1) // 2] += 1
    assert functions._mwu_null_counts(n1, n2) == counts
    assert counts == counts[::-1]  # symmetric, as the two-sided p-value assumes


def test_mannwhitneyu_matches_the_published_example():
    # SciPy's documented example: females against males, exact, "less": U1 = 3,
    # p = 0.05555555555555555 (= 7/126, the 7 splits with U2 >= 17 of C(9, 4)).
    males = [19.0, 22.0, 16.0, 29.0, 24.0]
    females = [20.0, 11.0, 17.0, 12.0]
    with mpmath.workdps(50):
        u, p = functions._mp_mannwhitneyu(mp(*females), mp(*males), mpmath.mpf(1), mpmath.mpf(-1))
        assert p == mpmath.mpf(7) / 126
    assert u == 3
    assert functions._mwu_exact_sf(2 * 17, 5, 4) == (7, 126)
    assert float(p) == 0.05555555555555555


def test_mannwhitneyu_ranks_and_statistic_match_counting():
    # U1 = #(x > y) + #(x = y) / 2 over all pairs, for every grid sample (ties included).
    for x, y, method, alt in functions._mwu_grid():
        if method != 0 or alt != 0:
            continue
        twice_pairs = sum(2 if a > b else 1 if a == b else 0 for a in x for b in y)
        assert functions._mwu_twice_u1(x, y) == twice_pairs


@pytest.mark.parametrize(
    "x, y",
    [((0, 0, 1, 2), (1, 1, 2, 3, 3)), ((5, 5, 5), (5, 6, 6, 7)), ((1, 2, 2, 2, 3, 4), (2, 4, 4))],
)
def test_mannwhitneyu_tie_corrected_variance_matches_enumeration(x, y):
    # s^2 against the variance of U1 over every split of the pooled values (average
    # ranks), which is the exact null variance with ties.
    xy = x + y
    n1, n = len(x), len(x) + len(y)
    us = [
        Fraction(
            functions._mwu_twice_u1([xy[i] for i in s], [xy[i] for i in range(n) if i not in s]), 2
        )
        for s in itertools.combinations(range(n), n1)
    ]
    mean = sum(us, Fraction(0)) / len(us)
    var = sum(((u - mean) ** 2 for u in us), Fraction(0)) / len(us)
    assert mean == Fraction(n1 * len(y), 2)
    with mpmath.workdps(50):
        # z without continuity correction at U = mean + 1 is 1 / s.
        twice_u = int(2 * (mean + 1))
        z = functions._mwu_z(twice_u, n1, len(y), functions._mwu_tie_term(x, y), False)
        assert close(1 / z**2, mpmath.mpf(var.numerator) / var.denominator)


@pytest.mark.parametrize("alt", [-1.0, 0.0, 1.0])
def test_mannwhitneyu_asymptotic_pvalue_matches_the_normal_cdf(alt):
    # The erfc tail against mpmath.ncdf at the tested U (U1, U2 or the larger), clipped at 1.
    x, y = functions._mwu_samples(10, 25, 150.5, False)
    u1 = Fraction(functions._mwu_twice_u1(x, y), 2)
    with mpmath.workdps(50):
        for method in (2.0, 3.0):
            u, p = functions._mp_mannwhitneyu(mp(*x), mp(*y), mpmath.mpf(method), mpmath.mpf(alt))
            assert u == mpmath.mpf(u1.numerator) / u1.denominator
            tested = {1.0: u1, -1.0: 250 - u1, 0.0: max(u1, 250 - u1)}[alt]
            mu, cc = mpmath.mpf(125), mpmath.mpf(0.5) if method == 2 else 0
            z = (mpmath.mpf(tested.numerator) / tested.denominator - mu - cc) / mpmath.sqrt(
                mpmath.mpf(10 * 25 * 36) / 12
            )
            expected = min((2 if alt == 0 else 1) * mpmath.ncdf(-z), 1)
            assert close(p, expected, digits=35)


def test_mannwhitneyu_automatic_method_follows_the_rule():
    # Exact when n1 <= 8 or n2 <= 8 and there are no ties; asymptotic with continuity
    # correction otherwise (SciPy's and statrs's rule).
    with mpmath.workdps(50):
        for x, y, method, alt in functions._mwu_grid():
            if method != 0:
                continue
            ties = len(set(x) | set(y)) < len(x) + len(y)
            rule = 2.0 if (len(x) > 8 and len(y) > 8) or ties else 1.0
            args = (mp(*x), mp(*y))
            auto = functions._mp_mannwhitneyu(*args, mpmath.mpf(0), mpmath.mpf(alt))
            assert auto == functions._mp_mannwhitneyu(*args, mpmath.mpf(rule), mpmath.mpf(alt))


def test_mannwhitneyu_rounding_bound_matches_the_tolerance_reason():
    with mpmath.workdps(30):
        worst, point = max(
            (functions.mannwhitneyu_rounding_bound(*p), p) for p in functions._mwu_grid()
        )
    assert 8.5e-14 < worst <= functions.MANNWHITNEYU.tolerance
    assert (len(point[0]), len(point[1]), point[2]) == (10, 25, 3.0)


def test_mannwhitneyu_two_sided_pvalue_is_capped_at_one():
    # U1 at its mean n1 n2 / 2: twice the tail is above 1 for both the exact and the
    # asymptotic methods (with continuity correction), so the p-value is exactly 1.
    x, y = (1.0, 4.0, 5.0, 8.0), (2.0, 3.0, 6.0, 7.0)  # U1 = 8 = 4 * 4 / 2
    with mpmath.workdps(50):
        for method in (1.0, 2.0):
            u, p = functions._mp_mannwhitneyu(mp(*x), mp(*y), mpmath.mpf(method), mpmath.mpf(0))
            assert u == 8
            assert p == 1


@pytest.mark.parametrize(
    "table", [(2, 1, 1, 2), (3, 1, 1, 3), (0, 3, 4, 1), (2, 3, 1, 4), (5, 0, 1, 2)]
)
def test_fisher_exact_counts_match_enumeration_of_labellings(table):
    # The hypergeometric counts against a different route: label every one of the
    # n1 + n2 items by row, choose every set of n of them for the first column, and
    # count the sets with x row-1 items. Every p-value is then a ratio of such counts.
    a, b, c, d = table
    n1, n = a + b, a + c
    rows = [0] * n1 + [1] * (c + d)
    counts = {}
    for chosen in itertools.combinations(range(len(rows)), n):
        x = sum(1 for i in chosen if rows[i] == 0)
        counts[x] = counts.get(x, 0) + 1
    w, lo, total = functions._fisher_weights(*table)
    assert {lo + i: v for i, v in enumerate(w)} == counts
    assert total == sum(counts.values())
    assert functions._fisher_pvalue_counts(*table, -1.0) == (
        sum(v for x, v in counts.items() if x <= a),
        total,
    )
    assert functions._fisher_pvalue_counts(*table, 1.0) == (
        sum(v for x, v in counts.items() if x >= a),
        total,
    )
    assert functions._fisher_pvalue_counts(*table, 0.0) == (
        sum(v for v in counts.values() if v <= counts[a]),
        total,
    )


def test_fisher_exact_greater_matches_less_on_the_second_column():
    # P(X >= a) is P(Y <= b) for Y the first-row count of the second column (the route
    # SciPy and statrs take), for every grid table.
    for a, b, c, d, alt in functions._fisher_grid():
        if alt != 1:
            continue
        a, b, c, d = int(a), int(b), int(c), int(d)
        greater = functions._fisher_pvalue_counts(a, b, c, d, 1.0)
        assert greater == functions._fisher_pvalue_counts(b, a, d, c, -1.0)


def test_fisher_exact_matches_the_published_examples():
    with mpmath.workdps(50):
        # SciPy's documentation: [[8, 2], [1, 5]], two-sided: odds ratio 20, p = 5/143.
        odds, p = functions._mp_fisher_exact(*mp(8, 2, 1, 5, 0))
        assert odds == 20 and p == mpmath.mpf(5) / 143
        # SciPy documents 0.034965034965034975, one ulp above 5/143 correctly rounded.
        assert float(p) == 0.03496503496503497
        assert math.nextafter(float(p), 1) == 0.034965034965034975
        # [[6, 2], [1, 4]]: x = 2, 6, 7 of the support 2..7, 0.0163 + 0.0816 + 0.00466.
        _, p = functions._mp_fisher_exact(*mp(6, 2, 1, 4, 0))
        assert p == mpmath.mpf(7 + 35 + 2) / 429
        assert abs(p - mpmath.mpf("0.10256")) < 1e-5
        # Fisher's tea tasting, 3 of 4 cups right: one-sided 17/70, two-sided 34/70.
        _, p = functions._mp_fisher_exact(*mp(3, 1, 1, 3, 1))
        assert p == mpmath.mpf(17) / 70
        _, p = functions._mp_fisher_exact(*mp(3, 1, 1, 3, 0))
        assert p == mpmath.mpf(34) / 70


def test_fisher_exact_pvalue_matches_mpmath_binomials():
    # The exact counts against hypergeometric probabilities from mpmath.binomial at
    # 60 digits, for the largest margins (701 terms).
    a, b, c, d = 380, 620, 320, 580
    with mpmath.workdps(60):
        n1, n2, n = a + b, c + d, a + c
        pmf = {
            x: mpmath.binomial(n1, x) * mpmath.binomial(n2, n - x) / mpmath.binomial(n1 + n2, n)
            for x in range(n + 1)
        }
        assert close(sum(pmf.values()), mpmath.mpf(1), digits=50)
        expected = {
            -1.0: sum(v for x, v in pmf.items() if x <= a),
            1.0: sum(v for x, v in pmf.items() if x >= a),
            0.0: sum(v for v in pmf.values() if v <= pmf[a] * (1 + mpmath.mpf(10) ** -40)),
        }
        for alt, e in expected.items():
            _, p = functions._mp_fisher_exact(*mp(a, b, c, d, alt))
            assert close(p, e, digits=45)


def test_fisher_exact_near_ties_leave_out_the_more_probable_table():
    # Each near-tie table has a table at most 1e-4 more probable, which the two-sided
    # p-value leaves out (it would be included under statrs's 1e-4 window).
    ties = functions._fisher_near_ties()
    grid = functions._fisher_grid()
    assert len(ties) == 12 and ties[0] == (6, 9, 16, 41)
    for table in ties:
        w, lo, _ = functions._fisher_weights(*table)
        wa = w[table[0] - lo]
        closer = [v for v in w if wa < v and (v - wa) * 10**4 <= v]
        assert closer
        count, _ = functions._fisher_pvalue_counts(*table, 0.0)
        assert count == sum(v for v in w if v <= wa)
        assert (*map(float, table), 0.0) in grid


def test_fisher_exact_rounding_bound_matches_the_tolerance_reason():
    worst, point = max(
        (functions.fisher_exact_rounding_bound(*p), p) for p in functions._fisher_grid()
    )
    assert 7.7e-14 < worst <= functions.FISHER_EXACT.tolerance
    a, b, c, d, _ = point
    assert (a + b, c + d, a + c) == (1000, 900, 700)
