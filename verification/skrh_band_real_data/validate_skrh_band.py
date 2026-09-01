"""Real-data verification: SKRH-BAND, the first pair this pipeline's FastMspec/single-taper
cross-spectrum computation AND the instrumented dispersion-curve picker have both been checked
against on real, large-scale data (N=10801, 1605 traces) -- not synthetic fixtures.

Requires three files not committed to this repo (bulk/derived real data, matching this project's
established practice -- see repo root README's "Getting the example SAC data" section for the
general pattern):
  - AFSKRH_XVBAND_win_3_all_matched_data.mat  (267MB; bluehive:
    /scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/test/
    processed_data/love/madagascar/)
  - SKRH_BAND_fastmspec.mat, SKRH_BAND_firstorder.mat  (164KB each; repovibranium:
    /volume1/web/FastMSPEC_data/madagascar_data/pre_computed_files/)
Pull via `ssh <host> "cat <path>"` (both hosts' SFTP/scp is restricted; plain `cat` over ssh
works around it), into this directory, before running.

Findings this script reproduces (full narrative: docs/notebook5_revamp_progress.md, Stage 3):
1. Cross-spectrum cross-validation against Sayan Swar's own precomputed MATLAB coh_sum:
   - FastMspec: ~3% relative L2 error, concentrated near coherence nulls -- the signature of the
     already-documented, intentionally-not-reproduced MATLAB complex-floor bug
     (python/ccf_pipeline/NOTES.md), now confirmed on real data for the first time.
   - single-taper: must apply detrend + 5% cosine taper before the plain-FFT coherency (this
     project's established "5% Cosine Single-Taper" technique) -- once applied, matches MATLAB to
     ~5e-9 relative error (machine precision). Skipping that preprocessing gives a large, easily
     misread ~46% mismatch that looks like a bug but isn't one.
2. Dispersion-curve picking, reproducing Sayan's own executed notebook
   (phasevel_compute_slide13and14.ipynb): FastMspec converges to a clean picked curve; single-taper
   does not, even after extra smoothing -- confirmed with the vendored+instrumented picker
   (python/dispcurve_pick), which additionally surfaces *why* (bad_quality_fraction,
   n_accepted_picks, freq_coverage_fraction) rather than a bare pass/fail.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.signal import savgol_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from dispcurve_pick import extract_dispcurve, DispersionCurveExceptionWithDiagnostics

HERE = Path(__file__).resolve().parent
REF_PATH = HERE.parents[1] / "data" / "reference" / "SDISPL.ASC"
DIST_KM = 290.295


def load_ref_curve():
    df = pd.read_csv(REF_PATH, sep=r"\s+", header=0)
    df.columns = [c.strip().replace("(", "_").replace(")", "").replace("/", "_") for c in df.columns]
    mode0 = df[df["LMODE"] == 0].sort_values("FREQUENCY_Hz")
    return mode0[["FREQUENCY_Hz", "C_KM_S"]].values


def coherence_from_mat(mat_path):
    """Sayan's own faxis-construction formula (phasevel_compute_slide13and14.ipynb cell 8),
    reproduced verbatim -- for odd T it actually builds a length-(T+1) array (a benign quirk in
    his own code: he always indexes both faxis and coh_sum with the same *integer* index array
    `ind = np.where(faxis > 0)[0]`, whose values all stay < T regardless). A boolean mask applied
    directly to coh_sum (length T, not T+1) breaks on that mismatch -- use integer indices, same
    as his own working code, not a boolean mask.
    """
    m = loadmat(mat_path)
    coh_sum, coh_num = m['coh_sum'].reshape(-1, 1), m['coh_num']
    T = len(coh_sum)
    faxis = np.concatenate([np.arange(0, (T - (T - 1) % 2) / 2 + 1), np.arange(-(T - T % 2) / 2, 0)]) / T
    ind = np.where(faxis > 0)[0]
    return faxis[ind].flatten(), np.real(coh_sum[ind] / coh_num).flatten()


def pick(label, faxis, coh, ref_curve, **kwargs):
    print(f"--- {label} ---")
    try:
        _, curve, diag = extract_dispcurve(
            faxis, coh, DIST_KM, ref_curve, cmin=1.2, cmax=4.8, filt_width=10, filt_height=1.0,
            x_step=0.05, pick_threshold=0, horizontal_polarization=False, manual_picking=False,
            plotting=False, return_diagnostics=True, **kwargs,
        )
        print(f"  CONVERGED: curve spans [{curve[0,0]:.4f}, {curve[-1,0]:.4f}] Hz, "
              f"velocity [{curve[:,1].min():.2f}, {curve[:,1].max():.2f}] km/s")
        print(f"  {diag}")
    except DispersionCurveExceptionWithDiagnostics as e:
        print(f"  DID NOT CONVERGE -- {e.diagnostics}")
    print()


def main():
    print("=== Part 1: cross-spectrum cross-validation vs. Sayan's own precomputed MATLAB coh_sum ===\n")
    raw = loadmat(HERE / "AFSKRH_XVBAND_win_3_all_matched_data.mat", simplify_cells=True)
    s1, s2 = raw['S1_data_mat'], raw['S2_data_mat']
    n_onesided = s1.shape[2] // 2 + 1

    r_fast = compute_crosscorr_mtc_fastmspec(s1, s2, wband=0.001, cutoff=1 - 1e-5, epsilon=1e-5)
    mat_fast = loadmat(HERE / "SKRH_BAND_fastmspec.mat")
    err_fast = np.linalg.norm(r_fast.coh_sum[:n_onesided] - mat_fast['coh_sum'].flatten()[:n_onesided]) / \
        np.linalg.norm(mat_fast['coh_sum'].flatten()[:n_onesided])
    print(f"FastMspec relative L2 error vs. MATLAB: {err_fast:.4e} "
          f"(expected ~3% -- known complex-floor bug, see ccf_pipeline/NOTES.md; not a regression)")

    s1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1))
    s2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2))
    fft_s1, fft_s2 = np.fft.fft(s1p, axis=2), np.fft.fft(s2p, axis=2)
    coh_trace = fft_s2 * np.conj(fft_s1)
    coh_trace = np.where(np.isnan(coh_trace / np.abs(fft_s1) / np.abs(fft_s2)), 0,
                          coh_trace / np.abs(fft_s1) / np.abs(fft_s2))
    coh_sum_single = coh_trace.sum(axis=(0, 1))
    mat_single = loadmat(HERE / "SKRH_BAND_firstorder.mat")
    err_single = np.linalg.norm(coh_sum_single[:n_onesided] - mat_single['coh_sum'].flatten()[:n_onesided]) / \
        np.linalg.norm(mat_single['coh_sum'].flatten()[:n_onesided])
    print(f"single-taper (detrend+5% cosine taper) relative L2 error vs. MATLAB: {err_single:.4e} "
          f"(expected ~5e-9 -- machine precision)\n")

    print("=== Part 2: dispersion-curve picking, reproducing Sayan's own SKRH-BAND result ===\n")
    ref_curve = load_ref_curve()
    faxis_fast, coh_fast = coherence_from_mat(HERE / "SKRH_BAND_fastmspec.mat")
    faxis_first, coh_first = coherence_from_mat(HERE / "SKRH_BAND_firstorder.mat")

    pick("FastMspec (raw)", faxis_fast, coh_fast, ref_curve, freqmin=0.01, freqmax=0.5)
    coh_first_smooth = savgol_filter(coh_first, window_length=20, polyorder=3)
    pick("Single-taper (smoothed, Sayan's own params)", faxis_first, coh_first_smooth, ref_curve,
         freqmin=0.05, freqmax=0.5)
    pick("Single-taper (raw, symmetric band -- fair comparison beyond Sayan's own notebook)",
         faxis_first, coh_first, ref_curve, freqmin=0.01, freqmax=0.5)


if __name__ == "__main__":
    main()
