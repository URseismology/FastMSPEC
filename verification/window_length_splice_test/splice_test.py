"""One-off test: does splicing several non-overlapping 3-hour windows into a genuinely longer,
continuous time series improve convergence for a representative Q4 (far) pair? Per direct
request, following the window-length recommendation in
docs/round2_hypothesis_evaluation.tex's Section "A concrete, quantitative fix".

Windowing convention confirmed from python/ccf_pipeline/prepare_data.py (ported from
legacy/matlab_source/lib/ccf_prepare_data_Z.m): windows have 50% overlap
(pts_begin = win_length*0.5*(iwin-1)+1+nstart), so window i and window i+2 are exactly
contiguous, non-overlapping. Array indices [0, 2, 4, ...] (0-based) within one day give a
genuine, non-duplicated splice.
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")

import numpy as np
import pandas as pd
from scipy.io import loadmat
from pathlib import Path

from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from dispcurve_pick import extract_dispcurve, load_reference_curve, build_template_family, \
    DispersionCurveExceptionWithDiagnostics

CUTOFF, EPSILON = 1 - 1e-5, 1e-5
PICK_FREQMIN, PICK_FREQMAX = 0.01, 0.5
PICK_CMIN, PICK_CMAX = 1.2, 4.8
CORRIDOR_KM_S, CORRIDOR_STEP_KM_S = 0.8, 0.05
TARGET_NW = 10.0  # matches quartile 1's present-day quality level, per the report's Table

MANIFEST = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/data/madagascar_stn_conn_ccflist.csv"
REF_CURVE = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/data/SDISPL.ASC"
# First attempt used only 10 days (coh_num=10 vs. the baseline's 1605) -- confounded the
# resolution test with a ~160x loss of multi-record averaging, and did not converge either way.
# Corrected: use ALL available days (matching the baseline's own full stack) so the only thing
# that changes between baseline and spliced is N, not coh_num.
N_DAYS_TO_SPLICE = None  # None means "use every available day"


def splice_day(s_day: np.ndarray) -> np.ndarray:
    """s_day: (n_window, n_samples) for one day. Returns one spliced 1-D trace, using windows
    at array indices [0, 2, 4, ...] (contiguous, non-overlapping per the 50%-overlap convention),
    dropping each window's first sample after the first (shared with the previous window's last
    sample) to avoid duplication.
    """
    idx = list(range(0, s_day.shape[0], 2))
    pieces = [s_day[idx[0]]] + [s_day[i][1:] for i in idx[1:]]
    return np.concatenate(pieces)


def run_one(net1, stn1, net2, stn2, label):
    manifest = pd.read_csv(MANIFEST)
    row = manifest[(manifest.net1 == net1) & (manifest.stn1 == stn1) &
                    (manifest.net2 == net2) & (manifest.stn2 == stn2)].iloc[0]
    dist_km = float(row.stndist)

    raw = loadmat(row.filelocation, simplify_cells=True)
    s1, s2 = raw["S1_data_mat"], raw["S2_data_mat"]
    if s1.ndim == 2:
        s1, s2 = s1[None, :, :], s2[None, :, :]
    n_day = s1.shape[0]

    # --- Baseline: standard 3-hour window, full coh_num stack (matches production) ---
    # compute_crosscorr_mtc_fastmspec requires exactly 3D (n_day, n_window, n_samples) input
    # (see _reshape_to_traces_by_samples) -- pass s1/s2 as-is, no manual reshape.
    n_samples_std = s1.shape[2]
    wband_std = TARGET_NW / n_samples_std
    r_std = compute_crosscorr_mtc_fastmspec(s1, s2, wband=wband_std, cutoff=CUTOFF, epsilon=EPSILON)
    faxis_std = np.fft.fftfreq(n_samples_std, d=1.0)
    pos_std = faxis_std > 0
    coh_pos_std = r_std.coh_sum[pos_std].real / r_std.coh_num
    print(f"{label}: baseline (3hr window, N={n_samples_std}, coh_num={r_std.coh_num}, "
          f"NW={TARGET_NW}, wband={wband_std:.6f})")
    result_std = pick(faxis_std[pos_std], coh_pos_std, dist_km, "baseline_3hr")

    # --- Spliced: N_DAYS_TO_SPLICE separate days, each spliced to one ~24hr continuous trace ---
    if N_DAYS_TO_SPLICE is None:
        day_indices = np.arange(n_day)
    else:
        day_indices = np.linspace(0, n_day - 1, N_DAYS_TO_SPLICE, dtype=int)
    spliced_traces_x, spliced_traces_y = [], []
    for d in day_indices:
        sx = splice_day(s1[d])
        sy = splice_day(s2[d])
        spliced_traces_x.append(sx)
        spliced_traces_y.append(sy)
    n_samples_spliced = len(spliced_traces_x[0])
    # (n_day, n_window, n_samples) shape required -- treat each spliced day as one "day" with
    # a single "window" (n_window=1), matching what compute_crosscorr_mtc_fastmspec expects.
    x_spliced = np.stack(spliced_traces_x, axis=0)[:, None, :]  # (n_traces, 1, n_samples_spliced)
    y_spliced = np.stack(spliced_traces_y, axis=0)[:, None, :]
    wband_spliced = TARGET_NW / n_samples_spliced
    r_spliced = compute_crosscorr_mtc_fastmspec(x_spliced, y_spliced, wband=wband_spliced,
                                                  cutoff=CUTOFF, epsilon=EPSILON)
    faxis_spliced = np.fft.fftfreq(n_samples_spliced, d=1.0)
    pos_spliced = faxis_spliced > 0
    coh_pos_spliced = r_spliced.coh_sum[pos_spliced].real / r_spliced.coh_num
    print(f"{label}: spliced ({N_DAYS_TO_SPLICE} days x ~24hr, N={n_samples_spliced}, "
          f"coh_num={r_spliced.coh_num}, NW={TARGET_NW}, wband={wband_spliced:.6f})")
    result_spliced = pick(faxis_spliced[pos_spliced], coh_pos_spliced, dist_km, "spliced_24hr")

    return result_std, result_spliced


def pick(faxis_pos, coh_pos, dist_km, tag):
    c_ref, f_lo, f_hi = load_reference_curve(Path(REF_CURVE), PICK_FREQMIN, PICK_FREQMAX)
    templates = build_template_family(c_ref, f_lo, f_hi, corridor_km_s=CORRIDOR_KM_S,
                                       step_km_s=CORRIDOR_STEP_KM_S)
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
        score = (diag.freq_coverage_fraction + 0.5 * (1 - min(diag.bad_quality_fraction, 1.0)) +
                 0.5 * min(diag.mean_amp_ratio / 5.0, 1.0)) if diag.converged else 0.0
        if score > best_score:
            best_score, best_delta, best_diag = score, delta, diag
    converged = bool(best_diag and best_diag.converged)
    print(f"  [{tag}] converged={converged} n_templates_converged={n_converged}/{len(templates)} "
          f"best_delta={best_delta} "
          f"bad_q={best_diag.bad_quality_fraction if best_diag else None} "
          f"cov={best_diag.freq_coverage_fraction if best_diag else None}")
    return {"converged": converged, "n_templates_converged": n_converged}


if __name__ == "__main__":
    run_one("AF", "SKRH", "XV", "BAEL", "AFSKRH_XVBAEL (Q4, 1013.7km, mean example)")
