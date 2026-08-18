"""Translation of firstNlambdaDPSS.m.

Computes the concentration eigenvalues of the first `first_n` DPSS tapers via
the same inverse-iteration approach as transitionDPSS.m. Only the eigenvalues
are returned, matching the original (the eigenvectors are computed internally
but discarded, same as the .m source).

Note: as in the original, the `k(1) > 1` branch that would extend the search
range is unreachable given the fixed starting range k = [1, first_n], and is
omitted here as dead code (see NOTES.md).

Two bug fixes vs. the .m source (both verified against scipy.signal.windows.dpss
to ~1e-15; see NOTES.md):

1. The original omits `v = v(end:-1:1);` before using v(j) as the
   inverse-iteration shift, unlike the equivalent block in transitionDPSS.m.
   Without the reversal, v(j) pairs the wrong tridiagonal eigenvalue with
   each seed vector.
2. The original's symmetric/antisymmetric polarization check uses
   `mod(k(2)+1-j,2)` -- the *descending*-index convention from
   transitionDPSS.m's upper-region loop -- while its seed vector uses the
   *ascending* `(j+k(1)-1)*t`, from that file's lower-region loop. These two
   conventions only agree by coincidence; mismatched, they cause the wrong
   symmetrize/antisymmetrize branch to be taken, which can badly corrupt an
   otherwise-converged eigenvector. Fixed by using the same index (seed_index)
   for both the seed and the parity check, consistent with transitionDPSS.m's
   lower-region loop that this function otherwise mirrors.
"""
from __future__ import annotations

import numpy as np

from ._dpss_transition import _diagonals
from ._tridiagonal import tridieig, tridisolve
from ._utils import fftfilt_autocorr


def first_n_lambda_dpss(n: int, w: float, first_n: int) -> np.ndarray:
    d, ee, t, s = _diagonals(n, w)
    k1, k2 = 1, first_n

    v = tridieig(d, np.concatenate([[0.0], ee]), n - k2 + 1, n - k1 + 1)
    v = v[::-1]  # bug fix vs. .m source, see module docstring

    lambda_all = np.zeros(k2 - k1 + 1)
    for j in range(1, k2 - k1 + 2):
        seed_index = j + k1 - 1
        e = np.sin(seed_index * t)
        e = tridisolve(ee, d - v[j - 1], e)
        e = tridisolve(ee, d - v[j - 1], e / np.linalg.norm(e))
        e = tridisolve(ee, d - v[j - 1], e / np.linalg.norm(e))

        parity_index = seed_index  # bug fix vs. .m source, see module docstring
        if parity_index % 2 == 0:
            if e[1] > 0:
                e = e - e[::-1]
            else:
                e = e[::-1] - e
        else:
            if e.sum() > 0:
                e = e + e[::-1]
            else:
                e = -e - e[::-1]

        e = e / np.linalg.norm(e)
        lambda_all[j - 1] = fftfilt_autocorr(e, s)

    return lambda_all
