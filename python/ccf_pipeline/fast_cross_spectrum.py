"""Batched (auto- or cross-) spectrum estimation using a fitted FastMultitaper,
translating mspec_fast.m's avgspec_sayan / avgspec_xy_sayan subfunctions
(functions/jSpectral/mspec_fast.m). These two MATLAB functions are
mathematically identical -- avgspec_sayan is just avgspec_xy_sayan called
with the same signal for both arguments -- so a single function here covers
both, matching what the MATLAB source does semantically even though it
keeps them as separate-named subfunctions.

This generalizes thomson_multitaper.FastMultitaper.spectral_estimate (single
signal, single trace, two-sided output) to: two signals (auto or cross),
many traces at once (batched), and MATLAB's one-sided output convention
(floor(N/2)+1 bins) -- because the caller (crosscorr_mtc.py) does its own
one-sided-to-two-sided reflection afterward, matching
ccf_compute_crosscorr_mtc_Z.m's structure.
"""
from __future__ import annotations

import numpy as np

from thomson_multitaper import FastMultitaper


def _complex_floor(z: np.ndarray, eps: float) -> np.ndarray:
    """Replicates MATLAB's `max(z, eps*maxz)` for possibly-complex z: MATLAB's
    max on complex arrays compares by magnitude and returns the actual
    (complex) element. maxz = the element of z with the largest magnitude.
    """
    flat_idx = np.argmax(np.abs(z))
    maxz = z.reshape(-1)[flat_idx]
    threshold = eps * maxz
    return np.where(np.abs(z) >= np.abs(threshold), z, threshold)


def fast_spectrum_batch(fmtse: FastMultitaper, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """One-sided (cross-)spectrum estimate for many traces at once.

    Parameters
    ----------
    fmtse : a fitted FastMultitaper(N, W, cutoff, epsilon) instance
    x, y : (N, n_traces) arrays, N == fmtse.n. Pass x is y (same array) for
        an auto-spectrum.

    Returns
    -------
    (floor(N/2)+1, n_traces) array -- complex in general (real up to
    round-off when x is y).
    """
    n = fmtse.n
    if x.shape[0] != n or y.shape[0] != n:
        raise ValueError(f"x/y must have {n} rows (fmtse.n), got {x.shape[0]}/{y.shape[0]}")
    eps = np.finfo(float).eps

    m = n
    l = m * int(np.ceil(2 * n / m))  # = 2n for this call convention (M=N always)
    vec_sinc = np.concatenate(
        [[2 * fmtse.w], fmtse.vec_half_sinc, np.zeros(l - 2 * n + 1), fmtse.vec_half_sinc[::-1]]
    )

    fx = np.fft.fft(x, n=l, axis=0)
    fy = np.fft.fft(y, n=l, axis=0)
    fxy = fx * np.conj(fy)
    z0 = np.fft.ifft(np.fft.fft(fxy, n=l, axis=0) * vec_sinc[:, None], n=l, axis=0)
    if l > m:
        z0 = z0[0 :: l // m]
    z0 = z0[: n // 2 + 1]

    sx = fmtse.S[:, None, :] * x[:, :, None]  # (n, n_traces, r)
    sy = fmtse.S[:, None, :] * y[:, :, None]
    mmatx = np.fft.fft(sx, n=n, axis=0)[: n // 2 + 1]  # one-sided: (Fbins, n_traces, r)
    mmaty = np.fft.fft(sy, n=n, axis=0)[: n // 2 + 1]
    eigspec = mmatx * np.conj(mmaty)
    z1 = eigspec[:, :, fmtse.index_plus].sum(axis=2) - eigspec[:, :, ~fmtse.index_plus].sum(axis=2)

    z = (z0 + z1) / fmtse.K
    if x is y:
        # auto-spectrum: a real, non-negative power -- the floor guards
        # against tiny negative numerical artifacts from the fast
        # sinc-convolution approximation (matches FastMultitaper.m's own
        # floor, which is only ever applied to an auto-spectrum).
        return _complex_floor(z, eps)
    # cross-spectrum: genuinely complex, with no non-negativity constraint.
    # mspec_fast.m's avgspec_xy_sayan applies the same floor here too, but
    # MATLAB's max() on a complex array compares by magnitude, so it
    # silently overwrites any bin near a coherence null -- exactly where
    # the phase is physically meaningful -- with a real value carrying the
    # *global* max's phase. See NOTES.md's "Known upstream bug" section.
    return z
