"""Reproductions of figures from Karnik, Romberg & Davenport, "Thomson's
Multitaper Method Revisited," IEEE Trans. Info. Theory 68(7), 2022.

Figs 1-3 and Table I have no MATLAB driver anywhere in the codebase (confirmed
by a full-tree grep during planning) -- these are fresh implementations built
directly from the paper's own equations and stated parameters, using the
already-verified `thomson_multitaper` library. Figs 4-8 and Table II *do* have
MATLAB drivers (`Comparison_ARMA_largescale_*.m`, `FastMultitaper_SpeedTest.m`
in `legacy/matlab_source`'s sibling copy under `sayan-swar-translation/`) --
those are ported here, at reduced scale (see NOTEBOOK 1's own markdown cells
for why: the paper's exact N=2^18 / 1000-trial / K=1047 scale was benchmarked
and found to take ~5 minutes per FastMultitaper precompute call alone, making
literal reproduction infeasible in this environment).
"""
from __future__ import annotations

import time

import numpy as np

from thomson_multitaper import FastMultitaper, Multitaper, MultitaperAdaptive
from thomson_multitaper._utils import dpss


# ---------------------------------------------------------------------------
# Fig 1: Slepian basis eigenvalue clustering
# ---------------------------------------------------------------------------

def fig1_eigenvalues(n: int = 10000, w: float = 1.0 / 100, k: int = 1000) -> np.ndarray:
    """First `k` DPSS concentration eigenvalues for signal length `n`,
    bandwidth `w`. Paper: N=10000, W=1/100, k=1000 (Fig. 1)."""
    _, lam = dpss(n, n * w, k)
    return lam


# ---------------------------------------------------------------------------
# Fig 2: spectral windows psi(f) for several taper counts
# ---------------------------------------------------------------------------

def spectral_window(tapers: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    """psi(f) = (1/K) sum_k |sum_n s_k[n] exp(-j 2 pi f n)|^2, evaluated by
    direct (non-FFT) summation so it can be sampled at arbitrary frequencies.

    tapers : (N, K) array, columns are DPSS tapers.
    freqs  : (M,) array of frequencies (cycles/sample).
    """
    n = tapers.shape[0]
    k = tapers.shape[1]
    nvec = np.arange(n)
    # (M, N) phase matrix @ (N, K) tapers -> (M, K)
    phase = np.exp(-2j * np.pi * np.outer(freqs, nvec))
    dtft = phase @ tapers
    return (np.abs(dtft) ** 2).sum(axis=1) / k


def fig2_spectral_windows(n: int = 2000, w: float = 1.0 / 100, ks=(39, 36, 32, 29), n_freqs: int = 4001):
    """Returns (freqs, {K: psi(f)}) for the four taper counts in Fig. 2 / Table I."""
    kmax = max(ks)
    tapers, _ = dpss(n, n * w, kmax)
    freqs = np.linspace(-0.5, 0.5, n_freqs)
    out = {}
    for k in ks:
        out[k] = spectral_window(tapers[:, :k], freqs)
    return freqs, out


def table1_leakage(n: int = 2000, w: float = 1.0 / 100, ks=(39, 36, 32, 29)):
    """Table I: max psi(f) outside [-W,W], and integral of psi(f) inside
    [-W,W] (== Sigma_K^(1), the mean of (1 - lambda_k) over the K tapers)."""
    kmax = max(ks)
    tapers, lam_full = dpss(n, n * w, kmax)
    freqs_out = np.linspace(w, 0.5, 4001)  # W <= |f| <= 1/2, by symmetry use positive side
    rows = []
    for k in ks:
        psi_out = spectral_window(tapers[:, :k], freqs_out)
        max_leak = float(psi_out.max())
        sigma1 = float(np.mean(1 - lam_full[:k]))
        rows.append({"K": k, "max_psi_outside_W": max_leak, "sigma1_integral_inside_W": sigma1})
    return rows


# ---------------------------------------------------------------------------
# Fig 3: four-narrowband high-dynamic-range leakage demo
# ---------------------------------------------------------------------------

def fig3_signal(n: int = 2000, seed: int = 0, oversample: int = 64):
    """Generates x ~ CN(0, R) with R the Toeplitz covariance matrix implied
    by the paper's four-narrowband PSD (its Lemma 5: R = int S(f) e_f e_f^H
    df over continuous f, not a covariance diagonal in the length-n DFT
    basis): S(f) = 1e3 on [0.18,0.22], 1e9 on [0.28,0.32], 1e2 on
    [0.38,0.42], 1e1 on [0.78,0.82], 1e0 elsewhere. Returns (x, S) with S
    evaluated on the length-n DFT grid f = 0..N-1 / N.

    Approximates that continuous-frequency process by synthesizing on a much
    finer grid (m = n*oversample) via the same frequency-domain-shaping
    idiom as `nb3_helpers.nlnm_synthetic`, then truncating to the first n
    samples -- i.e. genuinely windowing a finite piece out of a longer
    realization, rather than assigning independent values directly to the
    n analysis-grid bins. That truncation is what makes real spectral
    leakage possible: a signal built and analyzed on the *same* n-point grid
    is, by construction, exactly diagonal in the periodogram's own DFT
    basis, so no periodogram -- regardless of the true PSD's dynamic range
    -- can ever show leakage against it. Verified against the paper's own
    reported qualitative result for this example (periodogram and K=39 both
    swamp the three weaker sources; K=29 resolves all four)."""
    m = n * oversample
    fm = np.arange(m) / m
    s_m = np.ones(m)
    s_m[(fm >= 0.18) & (fm <= 0.22)] = 1e3
    s_m[(fm >= 0.28) & (fm <= 0.32)] = 1e9
    s_m[(fm >= 0.38) & (fm <= 0.42)] = 1e2
    s_m[(fm >= 0.78) & (fm <= 0.82)] = 1e1
    rng = np.random.default_rng(seed)
    w = (rng.standard_normal(m) + 1j * rng.standard_normal(m)) / np.sqrt(2)
    x_long = np.fft.ifft(np.fft.fft(w) * np.sqrt(s_m))
    x = x_long[:n]

    f = np.arange(n) / n
    s = np.ones(n)
    s[(f >= 0.18) & (f <= 0.22)] = 1e3
    s[(f >= 0.28) & (f <= 0.32)] = 1e9
    s[(f >= 0.38) & (f <= 0.42)] = 1e2
    s[(f >= 0.78) & (f <= 0.82)] = 1e1
    return x, s


def fig3_estimates(x: np.ndarray, n: int, w: float = 1.0 / 100, k_full: int = 39, k_trimmed: int = 29):
    periodogram = np.abs(np.fft.fft(x)) ** 2 / n
    mt_full = Multitaper(n, w, k_full).spectral_estimate(x)
    mt_trimmed = Multitaper(n, w, k_trimmed).spectral_estimate(x)
    return periodogram, mt_full, mt_trimmed


# ---------------------------------------------------------------------------
# Figs 4-6 / Table II: 8-method ARMA comparison (scaled down from N=2^18)
# ---------------------------------------------------------------------------

def arma_filter_coeffs():
    """Same ARMA(p,q) pole/zero configuration as Comparison_ARMA_largescale_*.m."""
    z_poles = np.concatenate([
        0.98 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.03),
        0.90 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.09),
        0.95 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.12),
        0.97 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.16),
        0.95 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.30),
        0.95 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.33),
    ])
    z_zeros = np.concatenate([
        1.2 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.15),
        1.2 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.20),
        1.1 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.23),
        1.1 * np.exp(np.array([-1, 1]) * 1j * 2 * np.pi * 0.26),
    ])
    b = np.poly(z_zeros)  # zeros -> numerator
    a = np.poly(z_poles)  # poles -> denominator
    return b, a


def arma_true_psd(n: int, b: np.ndarray, a: np.ndarray) -> np.ndarray:
    q = len(b) - 1
    p = len(a) - 1
    nn = np.arange(n)[:, None]
    num = np.exp(-1j * 2 * np.pi * nn * np.arange(q, -1, -1) / n) @ b
    den = np.exp(-1j * 2 * np.pi * nn * np.arange(p, -1, -1) / n) @ a
    return np.abs(num / den) ** 2


def arma_realization(n: int, b: np.ndarray, a: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One length-n realization of the complex ARMA(p,q) process, via the
    same 'generate 2N samples, discard the first N (transient)' scheme as the
    .m source."""
    from scipy.signal import lfilter

    w = (rng.standard_normal(2 * n) + 1j * rng.standard_normal(2 * n)) / np.sqrt(2)
    x = lfilter(b, a, w)
    return x[n:]


class ArmaComparisonParams:
    """Scaled-down replacement for the paper's N=2^18, W_small=1.25e-4,
    W_large=2e-3 configuration. See NOTEBOOK 1 markdown for the timing
    benchmarks that motivated this scale: at the paper's exact parameters, a
    single FastMultitaper precompute call took ~5 minutes, making the
    original 1000-trial experiment take days. Here, N and W are both reduced
    so that K_small/K_large keep a similar (smaller-but-representative) ratio
    to the paper's 64:1047 (~16x); a single trial across all 8 methods runs
    in ~7s, making 100 trials (~12 min) a feasible one-off run.
    """

    def __init__(self, n: int = 4096, w_small: float = 0.002, w_large: float = 0.016,
                 cutoff: float = 1 - 1e-9, epsilon: float = 1e-9):
        self.n = n
        self.w_small = w_small
        self.w_large = w_large
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.k_small = int(np.floor(2 * n * w_small)) - 1
        self.k_large = int(np.floor(2 * n * w_large)) - 1


METHOD_LABELS = [
    "1. Periodogram",
    "2. DPSS Tapered Periodogram",
    "3. Exact Multitaper, small W, K=2NW-1",
    "4. Exact Multitaper, small W, K=2NW-O(log NW)",
    "5. Approximate (Fast) Multitaper, large W, K=2NW-1",
    "6. Approximate (Fast) Multitaper, large W, K=2NW-O(log NW)",
    "7. Adaptive Multitaper, small W",
    "8. Adaptive Multitaper, large W",
]


def run_8_methods(x: np.ndarray, params: ArmaComparisonParams, time_it: bool = False):
    """Runs all 8 spectral estimation methods on one realization `x`.
    Returns (estimates, precompute_times, run_times), each a dict keyed by
    method index 1..8. Mirrors Comparison_ARMA_largescale_1trial.m /
    FastMultitaper_SpeedTest.m's method definitions."""
    n = params.n
    est, tpre, trun = {}, {}, {}

    t0 = time.time(); tpre[1] = 0.0
    t0 = time.time()
    est[1] = np.abs(np.fft.fft(x)) ** 2 / n
    trun[1] = time.time() - t0

    t0 = time.time()
    s0, _ = dpss(n, 4, 1)
    tpre[2] = time.time() - t0
    t0 = time.time()
    est[2] = np.abs(np.fft.fft(s0[:, 0] * x)) ** 2
    trun[2] = time.time() - t0

    t0 = time.time()
    mt3 = Multitaper(n, params.w_small, params.k_small)
    tpre[3] = time.time() - t0
    t0 = time.time()
    est[3] = mt3.spectral_estimate(x)
    trun[3] = time.time() - t0

    t0 = time.time()
    mt4 = Multitaper(n, params.w_small, params.cutoff)
    tpre[4] = time.time() - t0
    t0 = time.time()
    est[4] = mt4.spectral_estimate(x)
    trun[4] = time.time() - t0

    t0 = time.time()
    fmt5 = FastMultitaper(n, params.w_large, params.k_large, params.epsilon)
    tpre[5] = time.time() - t0
    t0 = time.time()
    est[5] = fmt5.spectral_estimate(x)
    trun[5] = time.time() - t0

    t0 = time.time()
    fmt6 = FastMultitaper(n, params.w_large, params.cutoff, params.epsilon)
    tpre[6] = time.time() - t0
    t0 = time.time()
    est[6] = fmt6.spectral_estimate(x)
    trun[6] = time.time() - t0

    t0 = time.time()
    amt7 = MultitaperAdaptive(n, params.w_small, params.k_small)
    tpre[7] = time.time() - t0
    t0 = time.time()
    est[7] = amt7.spectral_estimate(x)[0]
    trun[7] = time.time() - t0

    t0 = time.time()
    amt8 = MultitaperAdaptive(n, params.w_large, params.k_large)
    tpre[8] = time.time() - t0
    t0 = time.time()
    est[8] = amt8.spectral_estimate(x)[0]
    trun[8] = time.time() - t0

    return est, tpre, trun


def log_deviation(est: np.ndarray, true_psd: np.ndarray) -> np.ndarray:
    return np.abs(10 * np.log10(np.maximum(est, 1e-300)) - 10 * np.log10(true_psd))


# ---------------------------------------------------------------------------
# Figs 7-8 / Table II timing columns: exact vs. fast multitaper speed test
# ---------------------------------------------------------------------------

class SpeedTestParams:
    """Scaled-down replacement for FastMultitaper_SpeedTest.m's vecN (up to
    2^20) / 100 trials / 9 epsilons. See NOTEBOOK 1 markdown for the
    empirical scaling measurements (precompute time grows roughly linearly
    in N at fixed K, and grows with K too) that motivated capping N at 2^15
    and cutting epsilon/trial counts, while keeping the paper's own
    W = 0.08*N^(-1/5) bandwidth formula and cutoff=1-1e-3 unchanged."""

    def __init__(self, n_max_pow2: int = 15, n_min_pow2: int = 8,
                 n_max_exact: int = 4096, cutoff: float = 1 - 1e-3,
                 epsilons=(1e-4, 1e-8, 1e-12), num_trials: int = 8):
        self.vec_n = [2 ** k for k in range(n_min_pow2, n_max_pow2 + 1)]
        self.n_max_exact = n_max_exact
        self.cutoff = cutoff
        self.epsilons = list(epsilons)
        self.num_trials = num_trials

    def bandwidth(self, n: int) -> float:
        return 0.08 * n ** (-1 / 5)


def speed_test_trial(params: SpeedTestParams, rng: np.random.Generator):
    """One trial: for each N in params.vec_n, times exact Multitaper (if
    N <= n_max_exact) and FastMultitaper at each epsilon. Returns a dict of
    result rows, matching the columns FastMultitaper_SpeedTest.m records."""
    rows = []
    for n in params.vec_n:
        w = params.bandwidth(n)
        x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        x = x / np.linalg.norm(x)

        if n <= params.n_max_exact:
            t0 = time.time()
            mt = Multitaper(n, w, params.cutoff)
            t_pre_exact = time.time() - t0
            t0 = time.time()
            y_exact = mt.spectral_estimate(x)
            t_run_exact = time.time() - t0
            k_exact = mt.K
        else:
            t_pre_exact = t_run_exact = k_exact = None
            y_exact = None

        for eps in params.epsilons:
            t0 = time.time()
            fmt = FastMultitaper(n, w, params.cutoff, eps)
            t_pre_fast = time.time() - t0
            t0 = time.time()
            y_fast = fmt.spectral_estimate(x)
            t_run_fast = time.time() - t0

            max_err = float(np.max(np.abs(y_exact - y_fast))) if y_exact is not None else None

            rows.append({
                "N": n, "epsilon": eps, "K_exact": k_exact, "K_fast": fmt.K, "r_fast": fmt.r,
                "precompute_exact_s": t_pre_exact, "run_exact_s": t_run_exact,
                "precompute_fast_s": t_pre_fast, "run_fast_s": t_run_fast,
                "max_abs_error": max_err,
            })
    return rows
