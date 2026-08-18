"""Translation of MultitaperAdaptive.m: exact multitaper spectral estimator
with Thomson's adaptive (data-driven) taper weighting. Original Copyright
2019, Santhosh Karnik.

Deviations from the .m source (see NOTES.md):
  - The original's `cutoff > eps` branch calls `median(1, ceil(...)+10, N)`
    (missing the `[...]` that Multitaper.m uses around the same expression),
    which in MATLAB would error or misbehave rather than compute a median of
    three values. Translated here as the evident intent: `median([1, ..., N])`.
  - The original does not clamp K to be at least 1 (Multitaper.m does, via
    `K = max(K,1)`); clamped here too for robustness, since an all-zero
    taper set would otherwise crash spectral_estimate.
"""
from __future__ import annotations

import numpy as np

from ._utils import datawrap, dpss


class MultitaperAdaptive:
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

    def spectral_estimate(self, x: np.ndarray, m: int | None = None):
        """Adaptively-weighted multitaper PSD estimate.

        Returns (Sf, Skf, wkf, num_iter): the combined estimate, the
        per-taper single-taper estimates, the final weights, and the
        iteration count -- matching the .m source's four outputs.
        """
        if m is None:
            m = self.n
        x = np.asarray(x)

        if self.n > m:
            sx = self.S * x[:, None]
            skf = np.stack([datawrap(sx[:, k], m) for k in range(self.K)], axis=1)
        else:
            skf = self.S * x[:, None]

        skf = np.fft.fft(skf, n=m, axis=0)
        skf = skf.real**2 + skf.imag**2

        sigma2 = float(np.vdot(x, x).real) / self.n
        oneminuslambdasigma2 = np.tile((1 - self.lambda_) * sigma2, (m, 1))

        wkf = np.tile(self.lambda_, (m, 1))
        sf = (skf[:, 0] + skf[:, 1]) / 2
        sf_old = np.zeros(m)
        num_iter = 0
        while np.sum(np.abs(sf - sf_old)) > 5e-4 * sigma2 / m:
            num_iter += 1
            sf_old = sf
            wkf = (sf[:, None] ** 2 * self.lambda_[None, :]) / (
                sf[:, None] * self.lambda_[None, :] + oneminuslambdasigma2
            ) ** 2
            sf = np.sum(wkf * skf, axis=1) / np.sum(wkf, axis=1)

        return sf, skf, wkf, num_iter
