import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec

WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5

# (a) Sayan's own matched_data.mat, run through our already-validated FastMspec
sayan = loadmat(
    "/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/test/"
    "processed_data/love/madagascar/AFSKRH_XVBAND_win_3_all_matched_data.mat",
    simplify_cells=True,
)
s1, s2 = sayan["S1_data_mat"], sayan["S2_data_mat"]
print("Sayan matched_data shape:", s1.shape)
n_samples = s1.shape[2]
r_sayan = compute_crosscorr_mtc_fastmspec(s1, s2, wband=WBAND, cutoff=CUTOFF, epsilon=EPSILON)
faxis = np.fft.fftfreq(n_samples, d=1.0)
pos = faxis > 0
coh_sayan = r_sayan.coh_sum[pos].real / r_sayan.coh_num
print(f"Sayan-sourced: coh_num={r_sayan.coh_num}, range {coh_sayan.min():.4f} to {coh_sayan.max():.4f}")

# (b) our gvib.h5-sourced reprocessing, already computed
gvib = np.load("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/report_figures/gvib_skrh_band_coh.npz")
freq_gvib, coh_gvib = gvib["freq"], gvib["coh"]
print(f"gvib-sourced: n_day={gvib['n_day']}, dist_km={gvib['dist_km']:.1f}")

# (c) the original precomputed MATLAB coh_sum
mlab = loadmat("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/SKRH_BAND_fastmspec.mat")
coh_matlab_full = mlab["coh_sum"].flatten()
coh_num_matlab = int(mlab["coh_num"].flatten()[0])
n_onesided = n_samples // 2 + 1
coh_matlab = coh_matlab_full[pos[:len(coh_matlab_full)]].real / coh_num_matlab if len(coh_matlab_full) == n_samples else None
if coh_matlab is None:
    # coh_sum may already be one-sided; handle both cases
    print("MATLAB coh_sum length:", len(coh_matlab_full), "vs n_samples:", n_samples)
    coh_matlab = coh_matlab_full.real / coh_num_matlab

# quantitative comparison, Sayan-sourced (our own FastMspec) vs MATLAB (the ORIGINAL Stage 3 check)
rel_err_sayan_vs_matlab = np.linalg.norm(coh_sayan - coh_matlab[:len(coh_sayan)]) / np.linalg.norm(coh_matlab[:len(coh_sayan)])
print(f"Sayan-sourced vs MATLAB relative L2 error: {rel_err_sayan_vs_matlab:.4f} (Stage 3's own ~3% baseline)")

# gvib vs Sayan-sourced -- interpolate onto common frequency grid for comparison
from scipy.interpolate import interp1d
common_mask = (faxis[pos] >= freq_gvib.min()) & (faxis[pos] <= freq_gvib.max())
interp_gvib = interp1d(freq_gvib, coh_gvib, bounds_error=False, fill_value=0)
coh_gvib_on_sayan_grid = interp_gvib(faxis[pos])

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(faxis[pos], coh_sayan, color='k', lw=0.9, label="Sayan-sourced (matched_data.mat), our FastMspec")
ax.plot(faxis[pos], coh_matlab[:len(coh_sayan)], color='tab:blue', lw=0.7, alpha=0.7, label="Original MATLAB coh_sum")
ax.plot(freq_gvib, coh_gvib, color='tab:red', lw=0.9, alpha=0.85, label="gvib.h5-sourced (this test), our FastMspec")
ax.set_xlim(0, 0.35)
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Coherence (real part)")
ax.legend(fontsize=9)
ax.set_title("AF.SKRH-XV.BAND: gvib.h5 reprocessing vs. known-good references")
fig.tight_layout()
fig.savefig("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/report_figures/gvib_vs_reference_skrh_band.png", dpi=150)
print("saved comparison plot")
