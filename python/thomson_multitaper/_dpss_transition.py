"""Translation of transitionDPSS.m and transitionDPSS_modif.m.

Computes DPSS (Slepian) eigenvectors/eigenvalues in the "transition region"
around the concentration cutoff, via inverse iteration on the tridiagonal
commuting matrix T (Slepian, Grunbaum), rather than a dense eigendecomposition.
Original: Copyright 2017, Santhosh Karnik, "Fast Slepian Transform" toolbox.

Translation notes (see NOTES.md for the full list):
  - MATLAB is 1-based; loop/index bookkeeping below keeps the *value* of each
    index variable identical to the .m source and only subtracts 1 at the
    point of array access, to minimize transcription risk in this bisection-
    /inverse-iteration-heavy code.
  - `lastwarn`-based convergence-warning capture (an unused `msg` variable in
    the original) is dropped; it never affected control flow in the source.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from ._tridiagonal import tridieig, tridisolve
from ._utils import fftfilt_autocorr


@dataclass
class LambdaAll:
    lambda0: np.ndarray
    lambda1: np.ndarray


def _diagonals(n: int, w: float):
    idx = np.arange(n)  # 0-based, i.e. MATLAB's (0:N-1)
    d = ((n - 1 - 2 * idx) ** 2) * 0.25 * np.cos(2 * np.pi * w)
    i1 = np.arange(1, n)  # MATLAB's (1:N-1), length N-1
    ee = i1 * (n - i1) / 2  # MATLAB: (1:N-1).*(N-1:-1:1)/2
    t = idx / (n - 1) * np.pi
    s = np.concatenate([4 * w * np.sinc(2 * w * np.arange(n - 1, 0, -1)), [2 * w]])
    return d, ee, t, s


def _compute_one(j_1based: int, k1_1based: int, v_j: float, d, ee, t, s, descending_from=None):
    """One inverse-iteration eigenvector/eigenvalue pair.

    `descending_from`, if given, replaces `sin((j+k1-1)*t)` with
    `sin((descending_from - j + 1)*t)`, matching the second (upper-region)
    loop in the .m source which counts j upward but indexes the seed
    frequency downward from k(2).
    """
    n = len(d)
    if descending_from is None:
        seed_index = j_1based + k1_1based - 1
    else:
        seed_index = descending_from + 1 - j_1based
    e = np.sin(seed_index * t)
    e = tridisolve(ee, d - v_j, e)
    e = tridisolve(ee, d - v_j, e / np.linalg.norm(e))
    e = tridisolve(ee, d - v_j, e / np.linalg.norm(e))

    if seed_index % 2 == 0:
        # anti-symmetric DPSS
        if e[1] > 0:
            e = e - e[::-1]
        else:
            e = e[::-1] - e
    else:
        # symmetric DPSS
        if e.sum() > 0:
            e = e + e[::-1]
        else:
            e = -e - e[::-1]

    e = e / np.linalg.norm(e)
    lam = fftfilt_autocorr(e, s)
    return e, lam


def transition_dpss(n: int, w: float, epsilon=None, cutoff=None):
    """Direct translation of transitionDPSS.m.

    Returns (S, lam, lambda_all, lowerindex, upperindex) with S as an (N, r)
    array of transition-region DPSS tapers and lam the matching eigenvalues.
    lowerindex/upperindex are 1-based DPSS indices, matching the .m source,
    so callers translating other 1-based .m code (e.g. FastMultitaper.m) can
    keep using them unmodified; convert to 0-based only where you index into
    Python arrays with them.
    """
    eps = np.finfo(float).eps
    if epsilon is None:
        epsilon = (4 * eps * np.sqrt(n),) * 2
    epsilon = np.atleast_1d(epsilon)
    if len(epsilon) == 1:
        epsilon0 = epsilon1 = float(epsilon[0])
    else:
        epsilon0, epsilon1 = float(epsilon[0]), float(epsilon[1])

    if cutoff is None:
        cutoff = int(np.ceil(2 * n * w - 0.5))

    if cutoff >= 1:
        K = int(np.median([1, np.floor(cutoff), n - 1]))
    elif cutoff > eps:
        K = int(np.ceil(2 * n * w - 0.5))
        epsilon0 = min(epsilon0, cutoff)
        epsilon1 = min(epsilon1, 1 - cutoff)
    else:
        raise ValueError("cutoff must be >= 1 or > eps")

    d, ee, t, s = _diagonals(n, w)

    # --- Lower region: epsilon0 < lambda < 1/2 ---
    halfgap0 = np.log(max(n * np.sqrt(np.sin(2 * np.pi * w)), 2)) * np.log(1.0 / epsilon0) / np.pi**2 + 3
    k1, k2 = K + 1, min(n, int(np.ceil(2 * n * w - 0.5)) + int(np.ceil(halfgap0)))
    if k1 > k2:
        S0 = np.zeros((n, 0))
        lambda0 = np.zeros(0)
        index0 = 0
        flag = False
    else:
        S0 = np.zeros((n, k2 - k1 + 1))
        lambda0 = np.zeros(k2 - k1 + 1)
        idx = 0
        flag = True

    j = 0
    while flag:
        v = tridieig(d, np.concatenate([[0.0], ee]), n - k2 + 1, n - k1 + 1)
        v = v[::-1]

        for j in range(1, k2 - k1 + 2):
            e, lam = _compute_one(j, k1, v[j - 1], d, ee, t, s)
            S0[:, idx + j - 1] = e
            lambda0[idx + j - 1] = lam
            if lam <= epsilon0:
                flag = False
                index0 = idx + j - 1
                break

        if k2 == n:
            index0 = idx + j
            flag = False
        if k2 < n and flag:
            k1, k2 = k2 + 1, min(n, k2 + int(np.ceil(np.log(n))))
            idx = idx + j
            S0 = np.concatenate([S0, np.zeros((n, k2 - k1 + 1))], axis=1)
            lambda0 = np.concatenate([lambda0, np.zeros(k2 - k1 + 1)])

    # --- Upper region: 1/2 <= lambda < 1-epsilon1 ---
    halfgap1 = np.log(max(n * np.sqrt(np.sin(2 * np.pi * w)), 2)) * np.log(1.0 / epsilon1) / np.pi**2 + 3
    k1, k2 = max(1, int(np.ceil(2 * n * w + 0.5)) - int(np.ceil(halfgap1))), K
    if k1 > k2:
        S1 = np.zeros((n, 0))
        lambda1 = np.zeros(0)
        index1 = 0
        flag = False
    else:
        S1 = np.zeros((n, k2 - k1 + 1))
        lambda1 = np.zeros(k2 - k1 + 1)
        idx = 0
        flag = True

    j = 0
    while flag:
        v = tridieig(d, np.concatenate([[0.0], ee]), n - k2 + 1, n - k1 + 1)

        for j in range(1, k2 - k1 + 2):
            e, lam = _compute_one(j, k1, v[j - 1], d, ee, t, s, descending_from=k2)
            S1[:, idx + j - 1] = e
            lambda1[idx + j - 1] = lam
            if lam >= 1 - epsilon1:
                flag = False
                index1 = idx + j - 1
                break

        if k1 == 1:
            index1 = idx + j
            flag = False
        if k1 > 1 and flag:
            k1, k2 = max(1, k1 - int(np.ceil(np.log(n)))), k1 - 1
            idx = idx + j
            S1 = np.concatenate([S1, np.zeros((n, k2 - k1 + 1))], axis=1)
            lambda1 = np.concatenate([lambda1, np.zeros(k2 - k1 + 1)])

    lam = np.concatenate([lambda1[index1 - 1 :: -1], lambda0[:index0]])
    lambda_all = LambdaAll(lambda0=lambda0, lambda1=lambda1)
    S = np.concatenate([S1[:, index1 - 1 :: -1], S0[:, :index0]], axis=1)
    upperindex = K + index0
    lowerindex = K - index1 + 1

    return S, lam, lambda_all, lowerindex, upperindex


def transition_dpss_modif(n_old: int, w_old: float, epsilon=None, cutoff=None):
    """Translation of transitionDPSS_modif.m ("transitionDPSS_new" internally).

    CAUTION (see NOTES.md): the original computes DPSS at a fixed reference
    length N=512 and cubic-spline-interpolates eigenvectors back to n_old.
    The source contains an apparent bug -- `squared(e)` is not a MATLAB
    built-in and `S1b` is preallocated with 0 columns even when columns will
    be written into it -- both reproduced/fixed below as noted inline. This
    function has not been numerically validated against the original; treat
    it as a best-effort translation pending review.
    """
    eps = np.finfo(float).eps
    if epsilon is None:
        epsilon = (4 * eps * np.sqrt(n_old),) * 2
    epsilon = np.atleast_1d(epsilon)
    if len(epsilon) == 1:
        epsilon0 = epsilon1 = float(epsilon[0])
    else:
        epsilon0, epsilon1 = float(epsilon[0]), float(epsilon[1])

    if cutoff is None:
        cutoff = int(np.ceil(2 * n_old * w_old - 0.5))

    if cutoff >= 1:
        K = int(np.median([1, np.floor(cutoff), n_old - 1]))
    elif cutoff > eps:
        K = int(np.ceil(2 * n_old * w_old - 0.5))
        epsilon0 = min(epsilon0, cutoff)
        epsilon1 = min(epsilon1, 1 - cutoff)
    else:
        raise ValueError("cutoff must be >= 1 or > eps")

    p = n_old * w_old
    n = 512
    w = p / n

    x = np.arange(n) / (n - 1)
    xnew = np.arange(n_old) / (n_old - 1)

    d, ee, t, s = _diagonals(n, w)

    def interp_and_normalize(e_512):
        spline = CubicSpline(x, e_512)
        return spline(xnew)

    # --- Lower region ---
    halfgap0 = np.log(max(n * np.sqrt(np.sin(2 * np.pi * w)), 2)) * np.log(1.0 / epsilon0) / np.pi**2 + 3
    k1, k2 = K + 1, min(n, int(np.ceil(2 * n * w - 0.5)) + int(np.ceil(halfgap0)))
    if k1 > k2:
        S0b = np.zeros((n_old, 0))
        lambda0 = np.zeros(0)
        index0 = 0
        flag = False
    else:
        lambda0 = np.zeros(k2 - k1 + 1)
        S0b = np.zeros((n_old, k2 - k1 + 1))
        idx = 0
        flag = True

    j = 0
    while flag:
        v = tridieig(d, np.concatenate([[0.0], ee]), n - k2 + 1, n - k1 + 1)
        v = v[::-1]

        for j in range(1, k2 - k1 + 2):
            e, lam = _compute_one(j, k1, v[j - 1], d, ee, t, s)
            lambda0[idx + j - 1] = lam
            e_interp = interp_and_normalize(e)
            S0b[:, idx + j - 1] = e_interp / np.linalg.norm(e_interp)
            if lam <= epsilon0:
                flag = False
                index0 = idx + j - 1
                break

        if k2 == n:
            index0 = idx + j
            flag = False
        if k2 < n and flag:
            k1, k2 = k2 + 1, min(n, k2 + int(np.ceil(np.log(n))))
            idx = idx + j
            lambda0 = np.concatenate([lambda0, np.zeros(k2 - k1 + 1)])
            S0b = np.concatenate([S0b, np.zeros((n_old, k2 - k1 + 1))], axis=1)

    # --- Upper region ---
    halfgap1 = np.log(max(n * np.sqrt(np.sin(2 * np.pi * w)), 2)) * np.log(1.0 / epsilon1) / np.pi**2 + 3
    k1, k2 = max(1, int(np.ceil(2 * n * w + 0.5)) - int(np.ceil(halfgap1))), K
    if k1 > k2:
        S1b = np.zeros((n_old, 0))
        lambda1 = np.zeros(0)
        index1 = 0
        flag = False
    else:
        lambda1 = np.zeros(k2 - k1 + 1)
        # NOTE: original MATLAB preallocates S1b as zeros(Nold,0) here (a bug --
        # it should mirror S1's (k2-k1+1) column count like S0b above). Fixed here.
        S1b = np.zeros((n_old, k2 - k1 + 1))
        idx = 0
        flag = True

    j = 0
    while flag:
        v = tridieig(d, np.concatenate([[0.0], ee]), n - k2 + 1, n - k1 + 1)

        for j in range(1, k2 - k1 + 2):
            e, lam = _compute_one(j, k1, v[j - 1], d, ee, t, s, descending_from=k2)
            lambda1[idx + j - 1] = lam
            e_interp = interp_and_normalize(e)
            # NOTE: original uses `e/sqrt(sum(squared(e)))`; squared() is not a
            # MATLAB builtin, almost certainly meant `e.^2`, i.e. the L2 norm.
            S1b[:, idx + j - 1] = e_interp / np.linalg.norm(e_interp)
            if lam >= 1 - epsilon1:
                flag = False
                index1 = idx + j - 1
                break

        if k1 == 1:
            index1 = idx + j
            flag = False
        if k1 > 1 and flag:
            k1, k2 = max(1, k1 - int(np.ceil(np.log(n)))), k1 - 1
            idx = idx + j
            lambda1 = np.concatenate([lambda1, np.zeros(k2 - k1 + 1)])
            S1b = np.concatenate([S1b, np.zeros((n_old, k2 - k1 + 1))], axis=1)

    lam = np.concatenate([lambda1[index1 - 1 :: -1], lambda0[:index0]])
    lambda_all = LambdaAll(lambda0=lambda0, lambda1=lambda1)
    S = np.concatenate([S1b[:, index1 - 1 :: -1], S0b[:, :index0]], axis=1)
    upperindex = K + index0
    lowerindex = K - index1 + 1

    return S, lam, lambda_all, lowerindex, upperindex
