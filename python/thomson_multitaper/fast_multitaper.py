"""Translation of FastMultitaper.m: FFT-accelerated multitaper spectral
estimator using transition-region DPSS correction. Original Copyright 2019,
Santhosh Karnik.
"""
from __future__ import annotations

import numpy as np

from ._dpss_transition import transition_dpss
from ._utils import datawrap


class FastMultitaper:
    def __init__(self, n: int, w: float, cutoff: float, epsilon=None):
        self.n = n
        self.w = w
        self.cutoff = cutoff
        self.epsilon = epsilon
        eps = np.finfo(float).eps

        s, lam, lambda_all, lowerindex, upperindex = transition_dpss(n, w, epsilon, cutoff)

        if cutoff >= 1:
            k = cutoff
            eig_weights = np.concatenate([np.ones(k - lowerindex + 1), np.zeros(upperindex - k)]) - lam
        elif cutoff > eps:
            eig_weights = (lam > cutoff).astype(float) - lam
            k = int(np.count_nonzero(lam > cutoff)) + lowerindex - 1
        else:
            raise ValueError("cutoff must be >= 1 or > eps")

        self.S = s * np.sqrt(np.abs(eig_weights))[None, :]
        self.index_plus = eig_weights > 0
        self.r = len(eig_weights)
        self.K = k
        self.lambda_ = lam
        self.lambda_all = lambda_all

        self.vec_half_sinc = 2 * w * np.sinc(2 * w * np.arange(1, n))

    def spectral_estimate(self, x: np.ndarray, m: int | None = None) -> np.ndarray:
        if self.n != x.shape[0]:
            raise ValueError("x has incorrect length")
        if x.ndim > 1 and x.shape[1] > 1:
            raise ValueError("x must be a single column vector")

        if m is None:
            m = self.n
        eps = np.finfo(float).eps

        l = m * int(np.ceil(2 * self.n / m))
        vec_sinc = np.concatenate(
            [
                [2 * self.w],
                self.vec_half_sinc,
                np.zeros(l - 2 * self.n + 1),
                self.vec_half_sinc[::-1],
            ]
        )
        fx = np.fft.fft(x, n=l, axis=0)
        z0 = np.fft.ifft(np.fft.fft(fx.real**2 + fx.imag**2, n=l, axis=0) * vec_sinc, n=l, axis=0)
        if l > m:
            z0 = z0[0 :: l // m]

        if self.n > m:
            sx = self.S * x[:, None]
            a = np.stack([datawrap(sx[:, k], m) for k in range(self.r)], axis=1)
        else:
            a = self.S * x[:, None]

        a = np.fft.fft(a, n=m, axis=0)
        a = a * np.conj(a)
        z1 = a[:, self.index_plus].sum(axis=1) - a[:, ~self.index_plus].sum(axis=1)

        z = (z0 + z1) / self.K
        z = z.real  # z0/z1 are real up to FFT round-off; matches MATLAB's implicit real cast in max()
        max_z = z.max()
        return np.maximum(z, eps * max_z)
