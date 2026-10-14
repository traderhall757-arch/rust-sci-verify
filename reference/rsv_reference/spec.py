"""Description of one reference dataset (one mathematical function)."""

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import mpmath

# Tags are assigned only by rules in this package, never by hand.
POLE = "pole"  # the function has a pole or is undefined at this input
OVERFLOW = "overflow"  # the true value is finite but exceeds the largest double
UNSTABLE = "unstable-reference"  # mpmath gave different doubles at 50 and 100 digits

PRECISIONS = (50, 100)  # decimal digits; both must round to the same double


@dataclass(frozen=True)
class Spec:
    """A scalar function y = f(inputs) with its references and test grid."""

    name: str
    description: str
    inputs: tuple[str, ...]
    tolerance: float
    tolerance_reason: str
    grid_description: str
    grid: Callable[[], Sequence[tuple[float, ...]]]
    mpmath_fn: Callable[..., mpmath.mpf]
    scipy_name: str
    scipy_fn: Callable[..., float]
    # Input-based tag rules, e.g. "x is a non-positive integer" -> POLE.
    # Points tagged POLE are not evaluated with mpmath.
    input_tags: Callable[..., list[str]] = field(default=lambda *args: [])
