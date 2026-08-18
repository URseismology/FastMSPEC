"""Translation of ccf_compute_crosscorr_Z.m / _T.m -- the REAL entry point
of this pipeline (see docs/plan_ccf_mtc_translation.md): a dispatcher that
routes to one of three coherency computations depending on `filttype`,
called by the only two real driver scripts found in the codebase.

MATLAB dispatcher logic (identical in _Z.m and _T.m):
    if IsMultiTaper || IsFTN || IsOBN || IsSpecWhiten:
        # ccf_preprocess_filter_data.m already produced frequency-domain
        # data for one of these techniques; use it directly.
        fftS1, fftS2 = S1_data_mat, S2_data_mat
    elif IsMspec:
        return ccf_compute_crosscorr_mtc_Z(...)   # <- what this package translates
    else:
        fftS1, fftS2 = fft(S1_data_mat,3), fft(S2_data_mat,3)
    coh_trace = fftS2 .* conj(fftS1)  (note direction, Z convention)
    coh_trace = coh_trace ./ abs(fftS1) ./ abs(fftS2); nan->0
    -> ccf_save_computed_ccf_Z (multi-level day/month/full stacking)

Only the IsMspec path (this package's actual translation target, Phases
0-3) has been verified against Octave/real data. The other two branches
are implemented for completeness/fidelity to the dispatcher's structure,
but are NOT the exercised production path (every real config found sets
IsMspec=1 with the other three flags off) and have only been smoke-tested.
Multi-level stacking (ccf_save_computed_ccf_Z's day/month/single-stack
aggregation, as opposed to the full-stack-only save the _mtc_ functions
do internally) is not translated -- every real config found only uses
IsOutputFullstack=1 with the others off, and the _mtc_ path already
handles that case itself.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from thomson_multitaper import FastMultitaper

from .crosscorr_mtc import (
    CrosscorrResult,
    compute_crosscorr_mtc_fastmspec,
    compute_crosscorr_mtc_mspec,
    compute_crosscorr_mtc_mspecbestk,
)


@dataclass
class FilterConfig:
    """Mirrors the subset of MATLAB's `filttype` struct this dispatcher
    reads. Only one of is_multitaper/is_ftn/is_obn/is_specwhiten/is_mspec
    should be true, matching the .m source's implicit assumption (see
    ccf_preprocess_filter_data.m's NOTES.md caveat about these being
    mutually-exclusive "final transform" stages).
    """

    dt: float = 1.0
    is_multitaper: bool = False
    is_ftn: bool = False
    is_obn: bool = False
    is_specwhiten: bool = False
    is_mspec: bool = False
    technique: Literal["FastMspec", "Mspec", "MspecBestK"] | None = None
    wband: float | None = None
    cutoff: float | None = None
    epsilon: float | None = None
    nw: float | None = None
    k_taps: int | None = None
    fmtse: FastMultitaper | None = None


def _plain_coherency(fft_s1: np.ndarray, fft_s2: np.ndarray) -> tuple[np.ndarray, int]:
    """The non-IsMspec branches' shared tail: coh_trace = fftS2*conj(fftS1),
    normalized, NaN->0. Returns (coh_trace, coh_num) with coh_trace still
    3-D (day, window, freq) -- unlike the IsMspec path, which already sums
    over traces internally.
    """
    coh_trace = fft_s2 * np.conj(fft_s1)
    coh_num = coh_trace.shape[0] * coh_trace.shape[1]
    coh_trace = coh_trace / np.abs(fft_s1) / np.abs(fft_s2)
    coh_trace = np.where(np.isnan(coh_trace), 0, coh_trace)
    return coh_trace, coh_num


def compute_crosscorr(
    s1_data_mat: np.ndarray,
    s2_data_mat: np.ndarray,
    config: FilterConfig,
) -> CrosscorrResult | tuple[np.ndarray, int]:
    """Python port of ccf_compute_crosscorr_Z.m / _T.m's dispatcher.

    Returns a CrosscorrResult (summed full-stack coherency) for the
    IsMspec path, matching compute_crosscorr_mtc_*'s own return type; or
    a (coh_trace, coh_num) tuple for the other two branches, since those
    return an unsummed 3-D trace for potential multi-level stacking
    (not itself translated here -- see module docstring).
    """
    if config.is_multitaper or config.is_ftn or config.is_obn or config.is_specwhiten:
        fft_s1, fft_s2 = s1_data_mat, s2_data_mat
        return _plain_coherency(fft_s1, fft_s2)

    if config.is_mspec:
        if config.technique == "FastMspec":
            return compute_crosscorr_mtc_fastmspec(
                s1_data_mat, s2_data_mat, config.wband, config.cutoff, config.epsilon, fmtse=config.fmtse
            )
        if config.technique == "Mspec":
            return compute_crosscorr_mtc_mspec(
                s1_data_mat, s2_data_mat, config.wband, config.nw, config.k_taps, dt=config.dt
            )
        if config.technique == "MspecBestK":
            return compute_crosscorr_mtc_mspecbestk(
                s1_data_mat, s2_data_mat, config.wband, config.cutoff, config.epsilon, dt=config.dt
            )
        raise ValueError(f"Unknown technique: {config.technique!r}")

    fft_s1 = np.fft.fft(s1_data_mat, axis=2)
    fft_s2 = np.fft.fft(s2_data_mat, axis=2)
    return _plain_coherency(fft_s1, fft_s2)
