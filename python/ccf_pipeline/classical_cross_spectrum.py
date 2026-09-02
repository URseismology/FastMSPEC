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

# Chunk size over the taper (K) axis -- see the "Memory: chunked over K" note below.
# Picked so a single chunk's (N, n_traces, K_CHUNK) working set stays a few GB even at this
# project's largest real trace counts (~1600), not tuned further than that.
_DEFAULT_K_CHUNK = 8


def classical_spectrum_batch(x: np.ndarray, y: np.ndarray, tapers: np.ndarray, dt: float = 1.0,
                              k_chunk: int = _DEFAULT_K_CHUNK) -> np.ndarray:
    """One-sided (cross-)spectrum estimate, classical multitaper average.

    Parameters
    ----------
    x, y : (N, n_traces) arrays
    tapers : (N, K) array of DPSS tapers (e.g. from thomson_multitaper.dpss)
    dt : sample interval, matching the .m source's `avgspec(...) .* dt`
    k_chunk : tapers processed per batch (see "Memory: chunked over K" below);
        does not change the result (confirmed against the pre-chunking implementation
        to ~1e-15 relative error, the same floating-point-summation-order-only
        difference chunking always introduces), only peak memory.

    Returns
    -------
    (floor(N/2)+1, n_traces) array.

    Memory: chunked over K
    -----------------------
    The direct translation of `avgspec` (broadcast all K tapers against every trace at once,
    FFT, multiply, average) needs `(N, n_traces, K)`-shaped intermediates -- for Mspec's K=80 at
    this project's real Madagascar trace counts (~1600 traces), that is tens of GB *per array*,
    with several such arrays alive simultaneously (the pre-FFT real tapered traces, the post-FFT
    complex spectra for both x and y, and their product) -- confirmed to exceed even a generous
    64GB allocation on bluehive, killed by the OOM killer within ~2 minutes (i.e. failing at
    array construction, not partway through a long computation -- see
    docs/notebook5_revamp_progress.md's Stage 4 log for the concrete numbers). MspecBestK's much
    smaller K (~13-15) never hit this wall, which is why it went unnoticed until Mspec was
    actually run at real scale. Fixed by processing `k_chunk` tapers at a time and accumulating a
    running sum instead of materializing all K in memory together -- mathematically identical to
    `eigspec.mean(axis=2)` (that's exactly `eigspec.sum(axis=2) / K`, just accumulated
    incrementally), bounding peak memory to a `k_chunk`-sized working set regardless of K.
    """
    n = x.shape[0]
    n_traces = x.shape[1]
    k_total = tapers.shape[1]
    if y.shape[0] != n or tapers.shape[0] != n:
        raise ValueError("x, y, and tapers must all have N rows")

    n_onesided = n // 2 + 1
    eigspec_sum = np.zeros((n_onesided, n_traces), dtype=np.complex128)

    for k_start in range(0, k_total, k_chunk):
        tapers_chunk = tapers[:, k_start:k_start + k_chunk]  # (N, this_chunk_size)

        sx = tapers_chunk[:, None, :] * x[:, :, None]  # (N, n_traces, this_chunk_size)
        sy = tapers_chunk[:, None, :] * y[:, :, None]
        mmatx = np.fft.fft(sx, n=n, axis=0)[:n_onesided]
        mmaty = np.fft.fft(sy, n=n, axis=0)[:n_onesided]

        eigspec_sum += (mmatx * np.conj(mmaty)).sum(axis=2)

    return eigspec_sum / k_total * dt
