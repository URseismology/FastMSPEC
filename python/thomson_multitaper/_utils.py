"""Helpers standing in for MATLAB built-ins used throughout the original .m files."""
from __future__ import annotations

import numpy as np
from scipy.signal.windows import dpss as _scipy_dpss


def dpss(n: int, nw: float, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Discrete prolate spheroidal sequences, matching MATLAB's dpss(N, NW, K).

    Returns
    -------
    S : (n, k) array, columns are the tapers, most concentrated first
    lam : (k,) array, concentration eigenvalues (energy ratios), descending
    """
    tapers, ratios = _scipy_dpss(n, nw, Kmax=k, sym=True, norm=2, return_ratios=True)
    return tapers.T, ratios


def datawrap(x: np.ndarray, m: int) -> np.ndarray:
    """Time-domain alias x (length N > M) down to length M by wrapping, matching
    MATLAB's datawrap(x, M): y[i] = sum_{j: j % M == i} x[j].
    """
    x = np.asarray(x)
    n = x.shape[0]
    pad = (-n) % m
    if pad:
        x = np.concatenate([x, np.zeros(pad, dtype=x.dtype)])
    return x.reshape(-1, m).sum(axis=0)


def fftfilt_autocorr(taper: np.ndarray, sinc_kernel: np.ndarray) -> float:
    """Reproduces `fftfilt(flipud(e), e)' * s` from the .m source: the FIR
    filter output (impulse response = flipped taper, input = taper, causal,
    output truncated to len(taper)) dotted with the sinc kernel `s`.

    MATLAB's fftfilt(b, x) with len(b) == len(x) == N is equivalent to a full
    linear convolution truncated to the first N samples, i.e.
    np.convolve(b, x, mode='full')[:N].
    """
    n = len(taper)
    b = taper[::-1]
    y = np.convolve(b, taper, mode="full")[:n]
    return float(y @ sinc_kernel)
