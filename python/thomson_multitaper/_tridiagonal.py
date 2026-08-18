"""Tridiagonal linear algebra kernels.

Direct translations of tridisolve.m and tridieig.m from
ThomsonsMethodRevisitedExperiments/. The original tridisolve.m ships as a
MATLAB MEX-file stub with no visible source (its docstring cites Golub &
Van Loan, "Matrix Computations", 2nd ed., p.156); the implementation below
is the standard Thomas algorithm for a symmetric tridiagonal system, which
is what that reference describes. tridieig.m's source was available (MATLAB,
by C. Moler) and is translated near-literally, including its 1-based index
bookkeeping (kept via padding arrays) to minimize transcription risk in a
numerically delicate bisection routine.
"""
from __future__ import annotations

import numpy as np


def tridisolve(e: np.ndarray, d: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve A x = b where A = diag(e,-1) + diag(d,0) + diag(e,1).

    Parameters
    ----------
    e : off-diagonal, length N-1
    d : diagonal, length N
    b : right-hand side, length N

    Returns
    -------
    x : solution, length N

    Matches MATLAB's tridisolve(e, d, b, N) (the trailing N is redundant
    with len(d) and is dropped here). Assumes A is non-singular.
    """
    d = np.array(d, dtype=float, copy=True)
    x = np.array(b, dtype=float, copy=True)
    e = np.asarray(e, dtype=float)
    n = len(d)

    for k in range(1, n):
        mu = e[k - 1] / d[k - 1]
        d[k] = d[k] - mu * e[k - 1]
        x[k] = x[k] - mu * x[k - 1]

    x[n - 1] = x[n - 1] / d[n - 1]
    for k in range(n - 2, -1, -1):
        x[k] = (x[k] - e[k] * x[k + 1]) / d[k]

    return x


def tridieig(c: np.ndarray, b: np.ndarray, m1: int, m2: int, eps1: float = 0.0) -> np.ndarray:
    """Eigenvalues with (1-based) indices m1..m2 of a symmetric tridiagonal matrix.

    A = diag(b[1:], -1) + diag(c, 0) + diag(b[1:], 1); b[0] is ignored, matching
    the MATLAB convention `b(1) = 0`. Uses Sturm-sequence bisection (Wilkinson's
    algorithm), same as MATLAB's tridieig.m (C. Moler, MathWorks).

    Parameters
    ----------
    c : diagonal, length n
    b : off-diagonal with a leading placeholder, length n (b[0] unused)
    m1, m2 : 1-based inclusive index range of eigenvalues to compute, ascending
    eps1 : absolute tolerance; if <= 0, a machine-precision default is used

    Returns
    -------
    eigenvalues with indices m1..m2, ascending, as a 1-D array of length m2-m1+1
    """
    eps = np.finfo(float).eps
    n = len(c)

    # Internal arrays are 1-indexed via padding (index 0 unused) to mirror the
    # MATLAB source exactly; this is deliberate, not an oversight.
    c1 = np.zeros(n + 1)
    c1[1:] = c
    b1 = np.zeros(n + 1)
    b1[1:] = b
    b1[1] = 0.0
    beta = b1 * b1

    xmin = min(
        c1[n] - abs(b1[n]),
        min(c1[i] - abs(b1[i]) - (abs(b1[i + 1]) if i + 1 <= n else 0.0) for i in range(1, n)),
    )
    xmax = max(
        c1[n] + abs(b1[n]),
        max(c1[i] + abs(b1[i]) + (abs(b1[i + 1]) if i + 1 <= n else 0.0) for i in range(1, n)),
    )
    eps2 = eps * max(xmax, -xmin)
    if eps1 <= 0:
        eps1 = eps2
    eps2 = 0.5 * eps1 + 7 * eps2

    x0 = xmax
    x = np.zeros(n + 1)
    wu = np.zeros(n + 1)
    x[m1 : m2 + 1] = xmax
    wu[m1 : m2 + 1] = xmin

    for k in range(m2, m1 - 1, -1):
        xu = xmin
        for i in range(k, m1 - 1, -1):
            if xu < wu[i]:
                xu = wu[i]
                break
        if x0 > x[k]:
            x0 = x[k]
        while True:
            x1 = (xu + x0) / 2
            if x0 - xu <= 2 * eps * (abs(xu) + abs(x0)) + eps1:
                break
            a = 0
            q = 1.0
            for i in range(1, n + 1):
                if q != 0:
                    s = beta[i] / q
                else:
                    s = abs(b1[i]) / eps
                q = c1[i] - x1 - s
                if q < 0:
                    a += 1
            if a < k:
                if a < m1:
                    xu = x1
                    wu[m1] = x1
                else:
                    xu = x1
                    wu[a + 1] = x1
                    if x[a] > x1:
                        x[a] = x1
            else:
                x0 = x1
        x[k] = (x0 + xu) / 2

    return x[m1 : m2 + 1]
