"""Reference-value generator for rust-sci-verify.

mpmath, evaluated at two working precisions, is the arbiter; SciPy is recorded
alongside it as the practical reference. See docs/data-format.md for the file
format and docs/conventions.md for the classification rules.
"""

import os

# NumPy chooses SIMD code by CPU at import time. Its AVX-512 ("X86_V4") loops
# for exp and friends differ from the C library in the last bit, so SciPy
# functions built on them (e.g. scipy.stats.expon.sf) would give different
# reference files on different machines. Pin every machine to the AVX2 level
# that the CI runners use. This must run before NumPy is first imported;
# generate.check_numpy_dispatch() fails loudly if it did not take effect.
# Target names are those of the pinned NumPy (2.5).
NUMPY_DISABLED_TARGETS = "X86_V4 AVX512_ICL AVX512_SPR"
os.environ["NPY_DISABLE_CPU_FEATURES"] = NUMPY_DISABLED_TARGETS

# SciPy and NumPy call OpenBLAS (for example in scipy.stats.f_oneway), which
# also picks its kernels by CPU: its SkylakeX (AVX-512) kernels round
# differently from the Haswell (AVX2) kernels that the CI runners use, and the
# reference files must not depend on the machine. Pin every machine to Haswell,
# again before NumPy or SciPy load OpenBLAS; generate.check_numpy_dispatch()
# checks that the variable is in place.
OPENBLAS_CORETYPE = "Haswell"
os.environ["OPENBLAS_CORETYPE"] = OPENBLAS_CORETYPE
