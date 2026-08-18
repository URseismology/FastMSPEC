"""Translation of lib/ccf_compute_crosscorr_mtc_Z.m and _T.m.

Phase 1: the 'FastMspec' technique. Phase 2: 'Mspec' and 'MspecBestK'.
Phase 3: support for _T.m's calling-convention differences (see each
function's docstring) via optional parameters, rather than separate
duplicated _t functions -- the underlying math is identical between Z
and T; only how tapers are selected/injected differs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from thomson_multitaper import FastMultitaper, dpss

from .classical_cross_spectrum import classical_spectrum_batch
from .fast_cross_spectrum import fast_spectrum_batch


@dataclass
class CrosscorrResult:
    coh_sum: np.ndarray  # two-sided, length N (the window's sample count)
    coh_num: int  # number of traces (day*window) stacked
    taper_size: int  # number of tapers actually used


def _reshape_to_traces_by_samples(data_mat: np.ndarray) -> np.ndarray:
    """Reproduces MATLAB's `reshape(A, [], N1)` where A is (day, window,
    samples) and N1 = samples, followed by the transpose `'` used at the
    mspec_fast call site -- i.e. returns (samples, traces) with traces
    ordered day-fastest-then-window, matching MATLAB's column-major
    flattening. See docs/plan_ccf_mtc_translation.md's Phase 0 findings:
    this was empirically confirmed against Octave, not just reasoned about.
    """
    n_day, n_win, n_samp = data_mat.shape
    # (day, win, samp) -> (samp, win, day) -> reshape flattens last axis (day) fastest
    return data_mat.transpose(2, 1, 0).reshape(n_samp, n_day * n_win)


def _reflect_onesided_to_twosided(coh_sum_onesided: np.ndarray, n: int) -> np.ndarray:
    """Reconstructs the two-sided spectrum from the one-sided one, using the
    TRUE window length N's parity -- not the buggy `mod(length(coh_sum),2)`
    check in the original .m source, which uses the parity of the one-sided
    array's length instead of N. These agree only for N mod 4 in {1,2}; see
    docs/plan_ccf_mtc_translation.md and octave_verify_ccf/README.md for the
    full analysis (confirmed via Octave dry-run) of why the literal MATLAB
    check is wrong for N mod 4 in {0,3}.
    """
    if n % 2 == 0:
        neg = np.conj(np.flip(coh_sum_onesided[1:-1]))
    else:
        neg = np.conj(np.flip(coh_sum_onesided[1:]))
    return np.concatenate([coh_sum_onesided, neg])


def _coherency_stack(sxx: np.ndarray, syy: np.ndarray, sxy: np.ndarray, n_samples: int) -> np.ndarray:
    """Shared tail: cross-coherency, NaN->0, stack over traces, reflect
    one-sided->two-sided. Identical in ccf_compute_crosscorr_mtc_Z.m and
    _T.m for all three techniques.
    """
    coh = sxy / (np.sqrt(sxx) * np.sqrt(syy))
    coh = np.where(np.isnan(coh), 0, coh)
    coh_sum_onesided = coh.sum(axis=1)
    return _reflect_onesided_to_twosided(coh_sum_onesided, n_samples), coh.shape[1]


def compute_crosscorr_mtc_fastmspec(
    s1_data_mat: np.ndarray,
    s2_data_mat: np.ndarray,
    wband: float | None = None,
    cutoff: float | None = None,
    epsilon=None,
    fmtse: FastMultitaper | None = None,
) -> CrosscorrResult:
    """Python port of the 'FastMspec' branch, shared by _Z.m and _T.m.

    Z-style usage: pass wband/cutoff/epsilon, FastMultitaper is built inline
    (matches ccf_compute_crosscorr_mtc_Z.m, which always constructs it
    itself).
    T-style usage: pass a pre-built `fmtse` instead (matches
    ccf_compute_crosscorr_mtc_T.m, which uses `filttype.FMTSE` -- an object
    built once by the caller and reused across station pairs for
    efficiency, per docs/plan_ccf_mtc_translation.md Phase 3).

    Note the S2-before-S1 argument order in the original MATLAB (only in
    the 'FastMspec' branch, both Z and T), replicated here: SXY is computed
    as spectrum(S2, S1), i.e. X=S2, Y=S1 in the underlying cross-spectrum
    convention. This affects the conjugation/phase sense of the result and
    must not be "cleaned up" without re-deriving the expected sign.
    """
    n_samples = s1_data_mat.shape[2]
    if s2_data_mat.shape[2] != n_samples:
        raise ValueError("s1_data_mat and s2_data_mat must have the same window length")

    s1 = _reshape_to_traces_by_samples(s1_data_mat)  # (samples, traces)
    s2 = _reshape_to_traces_by_samples(s2_data_mat)

    if fmtse is None:
        if wband is None or cutoff is None or epsilon is None:
            raise ValueError("Provide either fmtse, or wband/cutoff/epsilon to build one")
        fmtse = FastMultitaper(n_samples, wband, cutoff, epsilon)
    elif fmtse.n != n_samples:
        raise ValueError(f"fmtse was built for N={fmtse.n}, but windows have {n_samples} samples")

    sxx = fast_spectrum_batch(fmtse, s2, s2)
    syy = fast_spectrum_batch(fmtse, s1, s1)
    sxy = fast_spectrum_batch(fmtse, s2, s1)

    coh_sum, coh_num = _coherency_stack(sxx, syy, sxy, n_samples)
    return CrosscorrResult(coh_sum=coh_sum, coh_num=coh_num, taper_size=fmtse.r)


def compute_crosscorr_mtc_mspec(
    s1_data_mat: np.ndarray,
    s2_data_mat: np.ndarray,
    wband: float | None = None,
    nw: float | None = None,
    k_taps: int | None = None,
    dt: float = 1.0,
) -> CrosscorrResult:
    """Python port of the 'Mspec' branch, shared by _Z.m and _T.m: classical
    multitaper cross-spectrum with plain DPSS tapers (no FastMultitaper
    fusion). Standard S1-before-S2 argument order (no swap, unlike
    'FastMspec', in both Z and T).

    Z-style usage: pass `wband`; taper count K is derived as
    `ceil(2*wband*N)-1` and the time-bandwidth product as `wband*N`,
    matching ccf_compute_crosscorr_mtc_Z.m's `sleptap(N1, Wband*N1,
    ceil(2*Wband*N1)-1)`.
    T-style usage: pass `nw`/`k_taps` directly (config-driven, not derived
    from a bandwidth), matching ccf_compute_crosscorr_mtc_T.m's
    `sleptap(N1, filttype.NW_mspec, filttype.K_taps_mspec)`.
    """
    n_samples = s1_data_mat.shape[2]
    if s2_data_mat.shape[2] != n_samples:
        raise ValueError("s1_data_mat and s2_data_mat must have the same window length")

    if nw is None or k_taps is None:
        if wband is None:
            raise ValueError("Provide either wband, or nw and k_taps")
        nw = wband * n_samples
        k_taps = int(np.ceil(2 * wband * n_samples)) - 1

    s1 = _reshape_to_traces_by_samples(s1_data_mat)
    s2 = _reshape_to_traces_by_samples(s2_data_mat)

    tapers, _ = dpss(n_samples, nw, k_taps)

    sxx = classical_spectrum_batch(s1, s1, tapers, dt)
    syy = classical_spectrum_batch(s2, s2, tapers, dt)
    sxy = classical_spectrum_batch(s1, s2, tapers, dt)

    coh_sum, coh_num = _coherency_stack(sxx, syy, sxy, n_samples)
    return CrosscorrResult(coh_sum=coh_sum, coh_num=coh_num, taper_size=tapers.shape[1])


def compute_crosscorr_mtc_mspecbestk(
    s1_data_mat: np.ndarray,
    s2_data_mat: np.ndarray,
    wband: float,
    cutoff: float,
    epsilon,
    dt: float = 1.0,
) -> CrosscorrResult:
    """Python port of the 'MspecBestK' branch -- identical between _Z.m and
    _T.m (no calling-convention differences for this technique). FastMultitaper
    used only to pick a taper count K, then plain DPSS tapers of that count
    (classical averaging, not the sinc-fused method). Standard S1-before-S2
    argument order.

    Note: the original MATLAB branch has a real bug -- it calls mspec_fast
    with only 4 output args (no totalMB), then unconditionally reads
    saved_ccf_path.psi_memory_space = totalMB afterward, which errors at
    runtime (confirmed by running the actual .m file in Octave). The
    coherency computation and its save() happen before that crash, so this
    Python port implements the (working, verifiable) coherency computation
    and simply omits the memory-diagnostic field entirely, rather than
    reproducing the crash.
    """
    n_samples = s1_data_mat.shape[2]
    if s2_data_mat.shape[2] != n_samples:
        raise ValueError("s1_data_mat and s2_data_mat must have the same window length")

    s1 = _reshape_to_traces_by_samples(s1_data_mat)
    s2 = _reshape_to_traces_by_samples(s2_data_mat)

    fmtse = FastMultitaper(n_samples, wband, cutoff, epsilon)
    tapers, _ = dpss(n_samples, int(np.ceil(n_samples * wband)), fmtse.K)

    sxx = classical_spectrum_batch(s1, s1, tapers, dt)
    syy = classical_spectrum_batch(s2, s2, tapers, dt)
    sxy = classical_spectrum_batch(s1, s2, tapers, dt)

    coh_sum, coh_num = _coherency_stack(sxx, syy, sxy, n_samples)
    return CrosscorrResult(coh_sum=coh_sum, coh_num=coh_num, taper_size=tapers.shape[1])
