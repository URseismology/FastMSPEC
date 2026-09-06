"""In-memory verification: recover AF.SKRH / XV.BAND raw traces directly from ADAMA_gvib.h5 via
python/dispcurve_pick/gvib_loader.py, compute the FastMspec cross-spectrum (no intermediate .mat
file -- everything held in numpy arrays), and compare against Stage 3's already-validated
SKRH-BAND result (verification/skrh_band_real_data/).

This script is now a thin driver over gvib_loader.build_pair_matched_data -- the windowing,
rotation, and zero-day-exclusion logic lives there (and is documented there in full, including
the real bugs this test caught and the fixes that followed). Keeping this file thin avoids the
duplication that let earlier versions of this script and the library drift apart silently.
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")

import numpy as np

from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from dispcurve_pick.gvib_loader import build_pair_matched_data

GVIB_PATH = "/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5"
# Fetch: curl -sSL https://raw.githubusercontent.com/URseismology/ADAMA/main/DataFiles/ADAMA_stalist.csv
STALIST = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/ADAMA_stalist.csv"
WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5

STA1, STA2 = "AF.SKRH", "XV.BAND"

print(f"Loading {STA1}/{STA2} pair-matched data from {GVIB_PATH} via gvib_loader ...")
pmd = build_pair_matched_data(GVIB_PATH, STA1, STA2, STALIST)
print(f"dist={pmd.dist_km:.1f} km (Stage 3's known value: 290.3 km -- cross-check), "
      f"az(1->2)={pmd.az_1to2_deg:.1f} deg, az(2->1)={pmd.az_2to1_deg:.1f} deg")
print(f"Assembled S1_data_mat/S2_data_mat: shape {pmd.S1_data_mat.shape} "
      f"({len(pmd.days_used)} days used), in memory, no .mat file written")

n_samples = pmd.S1_data_mat.shape[2]
r = compute_crosscorr_mtc_fastmspec(pmd.S1_data_mat, pmd.S2_data_mat, wband=WBAND, cutoff=CUTOFF,
                                     epsilon=EPSILON)
faxis = np.fft.fftfreq(n_samples, d=1.0)
pos = faxis > 0
coh = r.coh_sum[pos].real / r.coh_num
print(f"coh_num={r.coh_num}, coherence real-part range: {coh.min():.4f} to {coh.max():.4f}")

np.savez("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/report_figures/gvib_skrh_band_coh.npz",
         freq=faxis[pos], coh=coh, dist_km=pmd.dist_km, n_day=len(pmd.days_used))
print("Saved coherence for comparison against Stage 3's known-good result.")
