"""Each mpmath route used as the arbiter, checked against an independent one.

Agreement between 50 and 100 digits does not catch a wrong formula, so every
dataset's mpmath function is compared here with a different method (an exact
identity, a finite sum, or a different mpmath function) where both are valid.
"""

import json
import math
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
