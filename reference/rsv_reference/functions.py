"""The registry of reference datasets. One entry per mathematical function."""

import mpmath
import scipy.special as sc

from .spec import Spec


def _erfcx_grid() -> list[tuple[float]]:
    # Dense uniform part, step 0.25 (exact in binary), covering the
    # transition from the series/continued-fraction regions.
    dense = [i / 4 for i in range(-120, 121)]
    # Decades from 1e-300 to 1e300 on both sides. float("1eK") is correctly
    # rounded, so the grid does not depend on the platform's pow().
    decades = [float(f"1e{k}") for k in range(-300, 301, 10)]
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

ALL: tuple[Spec, ...] = (ERFCX,)
