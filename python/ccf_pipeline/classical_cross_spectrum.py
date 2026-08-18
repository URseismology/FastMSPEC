"""Batched (auto- or cross-) spectrum estimation using plain DPSS tapers,
translating mspec_fast.m's classical `avgspec` subfunction (the FMTSE=[]
code path) -- used by the 'Mspec' and 'MspecBestK' techniques, as opposed
to fast_cross_spectrum.py's fast_spectrum_batch (the FMTSE-fused path used
by 'FastMspec').

avgspec(mmat1, mmat2, N) = mean(mmat1 .* conj(mmat2), N) .* dt -- a simple
taper-averaged cross-periodogram, no sinc-kernel smoothing/fusion.
"""
from __future__ import annotations

import numpy as np


def classical_spectrum_batch(x: np.ndarray, y: np.ndarray, tapers: np.ndarray, dt: float = 1.0) -> np.ndarray:
    """One-sided (cross-)spectrum estimate, classical multitaper average.

    Parameters
    ----------
    x, y : (N, n_traces) arrays
    tapers : (N, K) array of DPSS tapers (e.g. from thomson_multitaper.dpss)
    dt : sample interval, matching the .m source's `avgspec(...) .* dt`

    Returns
    -------
    (floor(N/2)+1, n_traces) array.
    """
    n = x.shape[0]
    if y.shape[0] != n or tapers.shape[0] != n:
        raise ValueError("x, y, and tapers must all have N rows")

    sx = tapers[:, None, :] * x[:, :, None]  # (N, n_traces, K)
    sy = tapers[:, None, :] * y[:, :, None]
    mmatx = np.fft.fft(sx, n=n, axis=0)[: n // 2 + 1]
    mmaty = np.fft.fft(sy, n=n, axis=0)[: n // 2 + 1]

    eigspec = mmatx * np.conj(mmaty)
    return eigspec.mean(axis=2) * dt
