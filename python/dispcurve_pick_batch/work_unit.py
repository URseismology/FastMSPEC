"""Core per-(pair, technique) processing: load matched raw data, compute the cross-spectrum with
the right preprocessing for that technique, optionally cross-validate against Sayan Swar's own
precomputed MATLAB cross-spectrum, then scan the SDISPL.ASC +/- 0.8 km/s template family through
the instrumented picker and keep the best-scoring template as that work unit's pick.

Everything here was validated on real data first, at small scale, in
verification/skrh_band_real_data/ (Stage 3 of the Notebook 5 revamp) -- this module is that
validated logic, generalized from one hardcoded pair to any (Pair, technique) work unit. See
NOTES.md for the precomputed-MATLAB path pattern and why single-taper alone needs detrend+taper.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import (
    compute_crosscorr_mtc_fastmspec, compute_crosscorr_mtc_mspec, compute_crosscorr_mtc_mspecbestk,
)
from dispcurve_pick import (
    extract_dispcurve, DispersionCurveExceptionWithDiagnostics,
    load_reference_curve, build_template_family,
)

from .manifest import Pair

# Production config, matching legacy/matlab_source/entry_points/a2_ccf_run_crosscorr_T_mdg.m
WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5
NW_MSPEC, K_MSPEC = 100, 80

# Precomputed MATLAB cross-spectra live only for these two techniques, at this path pattern --
# confirmed by direct listing (Stage 4 investigation): 363 files each under a per-STA1 subdir.
MATLAB_TECHNIQUE_DIRS = {"FastMspec": "fastmspec", "single-taper": "firstorder"}
MATLAB_RESULTS_BASE = (
    "/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/"
    "results/test/love/madagascar/{technique_dir}/ccf/window3hr/fullStack/ccfTT/{sta1}/{sta1}_{sta2}_f.mat"
)

PICK_FREQMIN, PICK_FREQMAX = 0.01, 0.5
PICK_CMIN, PICK_CMAX = 1.2, 4.8
CORRIDOR_KM_S, CORRIDOR_STEP_KM_S = 0.8, 0.05


@dataclass
class WorkUnitResult:
    net1: str
    stn1: str
    net2: str
    stn2: str
    dist_km: float
    technique: str
    converged: bool
    best_delta_km_s: float | None
    bad_quality_fraction: float | None
    n_candidate_crossings: int | None
    n_accepted_picks: int | None
    freq_coverage_fraction: float | None
    mean_amp_ratio: float | None
    n_templates_converged: int
    n_templates_scanned: int
    matlab_coh_num: int | None
    matlab_rel_l2_error: float | None
    runtime_s: float
    error: str | None

    def as_dict(self) -> dict:
        return asdict(self)


def _coherency(s1_data_mat: np.ndarray, s2_data_mat: np.ndarray, technique: str):
    """Returns (coh_sum, coh_num) for the given technique. single-taper needs detrend + 5% cosine
    taper BEFORE the plain-FFT coherency (this project's established "5% Cosine Single-Taper"
    technique -- confirmed against real MATLAB output to 4.8e-9 relative error in Stage 3; skipping
    this preprocessing gives a misleading ~46% mismatch that looks like a bug but isn't one). The
    other three techniques use the un-preprocessed matched data directly, matching the production
    driver's IsDetrend=0/IsTaper=0 config for the IsMspec path.
    """
    if technique == "single-taper":
        s1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1_data_mat))
        s2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2_data_mat))
        fft_s1, fft_s2 = np.fft.fft(s1p, axis=2), np.fft.fft(s2p, axis=2)
        coh_trace = fft_s2 * np.conj(fft_s1)
        coh_num = coh_trace.shape[0] * coh_trace.shape[1]
        coh_trace = coh_trace / np.abs(fft_s1) / np.abs(fft_s2)
        coh_trace = np.where(np.isnan(coh_trace), 0, coh_trace)
        return coh_trace.sum(axis=(0, 1)), coh_num
    if technique == "FastMspec":
        r = compute_crosscorr_mtc_fastmspec(s1_data_mat, s2_data_mat, wband=WBAND, cutoff=CUTOFF, epsilon=EPSILON)
    elif technique == "Mspec":
        r = compute_crosscorr_mtc_mspec(s1_data_mat, s2_data_mat, nw=NW_MSPEC, k_taps=K_MSPEC, dt=1.0)
    elif technique == "MspecBestK":
        r = compute_crosscorr_mtc_mspecbestk(s1_data_mat, s2_data_mat, wband=WBAND, cutoff=CUTOFF, epsilon=EPSILON, dt=1.0)
    else:
        raise ValueError(f"Unknown technique: {technique!r}")
    return r.coh_sum, r.coh_num


def _matlab_cross_validate(pair: Pair, technique: str, coh_sum: np.ndarray, n_samples: int):
    """Returns (matlab_coh_num, matlab_rel_l2_error), or (None, None) if no precomputed MATLAB
    reference exists for this (pair, technique) -- only FastMspec and single-taper have one, and
    only for 363/380 pairs.
    """
    technique_dir = MATLAB_TECHNIQUE_DIRS.get(technique)
    if technique_dir is None:
        return None, None
    matlab_path = Path(MATLAB_RESULTS_BASE.format(technique_dir=technique_dir, sta1=pair.stn1, sta2=pair.stn2))
    if not matlab_path.exists():
        return None, None
    m = loadmat(matlab_path)
    coh_matlab = m["coh_sum"].flatten()
    coh_num_matlab = int(m["coh_num"].flatten()[0])
    n_onesided = n_samples // 2 + 1
    rel_err = float(
        np.linalg.norm(coh_sum[:n_onesided] - coh_matlab[:n_onesided]) /
        np.linalg.norm(coh_matlab[:n_onesided])
    )
    return coh_num_matlab, rel_err


def _score(diag) -> float:
    """Combined score for picking the best-fit template: convergence gates everything (a
    non-converged attempt scores 0 regardless of its partial diagnostics), then rewards frequency
    coverage, low bad-quality fraction, and a strong mean amplitude ratio -- the same 4 signals
    named in the revamp plan, combined as a simple weighted sum (not re-derived/tuned beyond that;
    a documented v1 choice, not a claim of optimality)."""
    if not diag.converged:
        return 0.0
    bad_q_term = 1.0 - min(diag.bad_quality_fraction, 1.0)
    amp_term = min(diag.mean_amp_ratio / 5.0, 1.0)  # saturate; ratios well above ~5 aren't meaningfully "more converged"
    return diag.freq_coverage_fraction + 0.5 * bad_q_term + 0.5 * amp_term


def process(pair: Pair, technique: str, ref_curve_path: Path) -> WorkUnitResult:
    t0 = time.time()
    try:
        raw = loadmat(pair.matched_data_path, simplify_cells=True)
        s1, s2 = raw["S1_data_mat"], raw["S2_data_mat"]
        # Pairs with only a single day of overlapping data collapse to a 2D (n_window, n_samples)
        # array in the .mat file (scipy/MATLAB drop the singleton day axis) instead of the usual
        # 3D (n_day, n_window, n_samples) -- found live, on GFOMA_XVLONA/GFOMA_XVMMBE (both ~2.4MB
        # files, vs. ~250MB+ for the heavily-sampled pairs tested so far), which crashed with
        # IndexError at s1.shape[2] before this fix. Restore the day axis explicitly rather than
        # assume every pair has multiple days.
        if s1.ndim == 2:
            s1 = s1[None, :, :]
        if s2.ndim == 2:
            s2 = s2[None, :, :]
        n_samples = s1.shape[2]
        dist_km = float(raw.get("stapairsinfo", {}).get("r", pair.dist_km))

        coh_sum, coh_num = _coherency(s1, s2, technique)
        matlab_coh_num, matlab_rel_err = _matlab_cross_validate(pair, technique, coh_sum, n_samples)

        faxis = np.fft.fftfreq(n_samples, d=1.0)
        pos = faxis > 0
        faxis_pos, coh_pos = faxis[pos], coh_sum[pos].real / coh_num

        c_ref, f_lo, f_hi = load_reference_curve(ref_curve_path, PICK_FREQMIN, PICK_FREQMAX)
        templates = build_template_family(c_ref, f_lo, f_hi, corridor_km_s=CORRIDOR_KM_S, step_km_s=CORRIDOR_STEP_KM_S)

        best_score, best_delta, best_diag = -1.0, None, None
        n_converged = 0
        for delta, c_template in templates.items():
            ref_curve_arr = np.column_stack([np.linspace(f_lo, f_hi, 200), c_template(np.linspace(f_lo, f_hi, 200))])
            try:
                _, _, diag = extract_dispcurve(
                    faxis_pos, coh_pos, dist_km, ref_curve_arr,
                    freqmin=f_lo, freqmax=f_hi, cmin=PICK_CMIN, cmax=PICK_CMAX,
                    filt_width=10, filt_height=1.0, x_step=0.05, pick_threshold=0,
                    horizontal_polarization=False, manual_picking=False, plotting=False,
                    return_diagnostics=True,
                )
            except DispersionCurveExceptionWithDiagnostics as e:
                diag = e.diagnostics
            if diag.converged:
                n_converged += 1
            score = _score(diag)
            if score > best_score:
                best_score, best_delta, best_diag = score, delta, diag

        return WorkUnitResult(
            net1=pair.net1, stn1=pair.stn1, net2=pair.net2, stn2=pair.stn2, dist_km=dist_km,
            technique=technique, converged=bool(best_diag and best_diag.converged),
            best_delta_km_s=best_delta if (best_diag and best_diag.converged) else None,
            bad_quality_fraction=best_diag.bad_quality_fraction if best_diag else None,
            n_candidate_crossings=best_diag.n_candidate_crossings if best_diag else None,
            n_accepted_picks=best_diag.n_accepted_picks if best_diag else None,
            freq_coverage_fraction=best_diag.freq_coverage_fraction if best_diag else None,
            mean_amp_ratio=best_diag.mean_amp_ratio if best_diag else None,
            n_templates_converged=n_converged, n_templates_scanned=len(templates),
            matlab_coh_num=matlab_coh_num, matlab_rel_l2_error=matlab_rel_err,
            runtime_s=time.time() - t0, error=None,
        )
    except Exception as e:  # noqa: BLE001 -- a batch of 1520 must never let one work unit's
        # exception kill the whole task; record it in the manifest instead. Re-raised nowhere:
        # every caller (run_plain.py, run_multiprocessing.py) treats this as a normal result row.
        return WorkUnitResult(
            net1=pair.net1, stn1=pair.stn1, net2=pair.net2, stn2=pair.stn2, dist_km=pair.dist_km,
            technique=technique, converged=False, best_delta_km_s=None, bad_quality_fraction=None,
            n_candidate_crossings=None, n_accepted_picks=None, freq_coverage_fraction=None,
            mean_amp_ratio=None, n_templates_converged=0, n_templates_scanned=0,
            matlab_coh_num=None, matlab_rel_l2_error=None, runtime_s=time.time() - t0,
            error=f"{type(e).__name__}: {e}",
        )
