"""Generate the coherence + seislib-picker diagnostic plot (coherence spectrum, KDE density,
zero crossings, dispersion curve) for one representative pair per distance quartile, for the
Round 2 hypothesis-evaluation report. Run on bluehive (raw .mat data lives there).
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
from scipy.io import loadmat
from pathlib import Path

from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from dispcurve_pick import extract_dispcurve, load_reference_curve, build_template_family, \
    DispersionCurveExceptionWithDiagnostics

WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5  # overridden per example below
N_SAMPLES = 10801
PICK_FREQMIN, PICK_FREQMAX = 0.01, 0.5
PICK_CMIN, PICK_CMAX = 1.2, 4.8
CORRIDOR_KM_S, CORRIDOR_STEP_KM_S = 0.8, 0.05

MANIFEST = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/data/madagascar_stn_conn_ccflist.csv"
REF_CURVE = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/data/SDISPL.ASC"
OUT_DIR = Path("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/report_figures")
OUT_DIR.mkdir(exist_ok=True)

EXAMPLES = [
    # (quartile, net1, stn1, net2, stn2, target_nw)
    (1, "XV", "BITY", "XV", "MAGY", 22.21),
    (2, "XV", "BAEL", "XV", "ZAKA", 22.18),
    (3, "XV", "BAND", "XV", "MAJA", 15.36),
    (4, "XV", "DGOS", "XV", "MAGY", 7.61),
]

manifest = pd.read_csv(MANIFEST)
c_ref, f_lo, f_hi = load_reference_curve(Path(REF_CURVE), PICK_FREQMIN, PICK_FREQMAX)
templates = build_template_family(c_ref, f_lo, f_hi, corridor_km_s=CORRIDOR_KM_S, step_km_s=CORRIDOR_STEP_KM_S)

for quartile, net1, stn1, net2, stn2, target_nw in EXAMPLES:
    row = manifest[(manifest.net1 == net1) & (manifest.stn1 == stn1) &
                    (manifest.net2 == net2) & (manifest.stn2 == stn2)].iloc[0]
    dist_km = float(row.stndist)
    wband = target_nw / N_SAMPLES

    raw = loadmat(row.filelocation, simplify_cells=True)
    s1, s2 = raw["S1_data_mat"], raw["S2_data_mat"]
    if s1.ndim == 2:
        s1, s2 = s1[None, :, :], s2[None, :, :]
    n_samples = s1.shape[2]

    r = compute_crosscorr_mtc_fastmspec(s1, s2, wband=wband, cutoff=CUTOFF, epsilon=EPSILON)
    coh_sum, coh_num = r.coh_sum, r.coh_num
    faxis = np.fft.fftfreq(n_samples, d=1.0)
    pos = faxis > 0
    faxis_pos, coh_pos = faxis[pos], coh_sum[pos].real / coh_num

    # Scan templates, keep the best-scoring one (matching work_unit.py's own selection logic)
    best_score, best_delta, best_diag = -1.0, None, None
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
            score = diag.freq_coverage_fraction + 0.5 * (1 - min(diag.bad_quality_fraction, 1.0)) + \
                    0.5 * min(diag.mean_amp_ratio / 5.0, 1.0)
        else:
            score = 0.0
        if score > best_score:
            best_score, best_delta = score, delta

    # Re-run the winning template WITH plotting on, saved to disk
    c_template = templates[best_delta]
    ref_curve_arr = np.column_stack([np.linspace(f_lo, f_hi, 200), c_template(np.linspace(f_lo, f_hi, 200))])
    savefig = str(OUT_DIR / f"q{quartile}_{net1}{stn1}_{net2}{stn2}.png")
    try:
        extract_dispcurve(
            faxis_pos, coh_pos, dist_km, ref_curve_arr,
            freqmin=f_lo, freqmax=f_hi, cmin=PICK_CMIN, cmax=PICK_CMAX,
            filt_width=10, filt_height=1.0, x_step=0.05, pick_threshold=0,
            horizontal_polarization=False, manual_picking=False, plotting=True,
            savefig=savefig, sta1=f"{net1}{stn1}", sta2=f"{net2}{stn2}",
        )
    except Exception as e:
        print(f"Q{quartile} {net1}{stn1}_{net2}{stn2}: plotting raised {type(e).__name__}: {e} "
              f"(figure may still have been saved before the raise)")
    print(f"Q{quartile} {net1}{stn1}_{net2}{stn2}: dist={dist_km:.1f}km NW={target_nw} "
          f"best_delta={best_delta} -> {savefig}")
