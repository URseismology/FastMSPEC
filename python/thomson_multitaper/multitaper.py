"""Translation of Multitaper.m: the exact (eigenvalue-cutoff) multitaper
spectral estimator. Original Copyright 2019, Santhosh Karnik.
"""
from __future__ import annotations

import numpy as np

from ._utils import datawrap, dpss


class Multitaper:
    """Exact multitaper spectral density estimator.

    Parameters
    ----------
    n : signal length the estimator is built for
    w : half-bandwidth (as a fraction of the sampling rate)
    cutoff : if >= 1, use the first `floor(cutoff)` DPSS tapers; if in
        (eps, 1), use all tapers whose concentration eigenvalue exceeds
        `cutoff`.
    """

    def __init__(self, n: int, w: float, cutoff: float):
        self.n = n
        self.w = w
        self.cutoff = cutoff
        eps = np.finfo(float).eps

        if cutoff >= 1:
            k = int(np.floor(cutoff + eps))
            s, lam = dpss(n, n * w, k)
        elif cutoff > eps:
            if cutoff > 0.5:
                k_est = int(np.ceil(2 * n * w))
            else:
                k_est = int(np.median([1, np.ceil(2 * n * w + (2 / np.pi**2) * np.log(n) * np.log(1.0 / cutoff)) + 10, n]))
            s_full, lam_full = dpss(n, n * w, k_est)
            k = max(int(np.count_nonzero(lam_full > cutoff)), 1)
            s, lam = s_full[:, :k], lam_full[:k]
        else:
            raise ValueError("cutoff must be >= 1 or > eps")

        self.S = s
        self.lambda_ = lam
        self.K = k

    def spectral_estimate(self, x: np.ndarray, m: int | None = None) -> np.ndarray:
        """y[m] = (1/K) sum_k |FFT_M(taper_k * x)[m]|^2, the multitaper PSD estimate."""
        if m is None:
            m = self.n
        x = np.asarray(x)

        if self.n > m:
            sx = self.S * x[:, None]
            a = np.stack([datawrap(sx[:, k], m) for k in range(self.K)], axis=1)
        else:
            a = self.S * x[:, None]

        a = np.fft.fft(a, n=m, axis=0)
        a = a.real**2 + a.imag**2
        return a.sum(axis=1) / self.K
