"""Section 2 precompute for Notebook 3: XV.BITY-XV.MAGY (Q1, 223 km) -- the
strongest converged FastMspec example from the revamp's own 380-pair + NW-sweep
testing. single-taper vs FastMspec raw coherence, current production bandwidth
(NW=22.21, the 1.5x NW_high value the sweep found best for this pair).

Built from ADAMA_gvib.h5 (same reproducible raw source as Section 1), not the
precomputed .mat -- so both worked examples in the notebook trace back to one
documented source.

Output: nb3_sec2_bity_magy.npz -- freq, coh_fastmspec, coh_single, dist_km,
coh_num_fm, coh_num_st, nw. (Dispersion-curve picking + all plotting happen live
in the notebook, on these 1-D arrays -- cheap, and keeps the hybrid reference
curve as a single source of truth there.)
"""
import json, sys, time
import numpy as np

D = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"
sys.path.insert(0, D + "/python")
from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from dispcurve_pick.gvib_loader import build_pair_matched_data

GVIB = "/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5"
STALIST = D + "/ADAMA_stalist.csv"
OUT = D + "/report_figures/nb3_sec2_bity_magy.npz"
TARGET_NW = 22.21
CUTOFF, EPSILON = 1 - 1e-5, 1e-5

print(f"Building XV.BITY-XV.MAGY from {GVIB} ...", flush=True)
pmd = build_pair_matched_data(GVIB, "XV.BITY", "XV.MAGY", STALIST)
s1, s2 = pmd.S1_data_mat, pmd.S2_data_mat
n_samples = s1.shape[2]
print(f"  dist={pmd.dist_km:.1f} km, shape={s1.shape}, {len(pmd.days_used)} days", flush=True)

faxis = np.fft.fftfreq(n_samples, d=1.0)
pos = faxis > 0

t0 = time.time()
r = compute_crosscorr_mtc_fastmspec(s1, s2, wband=TARGET_NW / n_samples, cutoff=CUTOFF, epsilon=EPSILON)
coh_fm = r.coh_sum[pos].real / r.coh_num
print(f"  FastMspec {time.time()-t0:.0f}s  K={r.taper_size}  coh_num={r.coh_num}", flush=True)

s1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1))
s2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2))
fft_s1, fft_s2 = np.fft.fft(s1p, axis=2), np.fft.fft(s2p, axis=2)
ct = fft_s2 * np.conj(fft_s1) / np.abs(fft_s1) / np.abs(fft_s2)
ct = np.where(np.isnan(ct), 0, ct)
coh_num_st = ct.shape[0] * ct.shape[1]
coh_st = ct.sum(axis=(0, 1))[pos].real / coh_num_st

np.savez(OUT, freq=faxis[pos], coh_fastmspec=coh_fm, coh_single=coh_st,
         dist_km=pmd.dist_km, coh_num_fm=r.coh_num, coh_num_st=coh_num_st,
         nw=TARGET_NW, taper_size_fm=r.taper_size)
print("DONE ->", OUT, flush=True)
print(json.dumps(dict(dist_km=round(float(pmd.dist_km), 1), coh_num_fm=int(r.coh_num),
                      coh_num_st=int(coh_num_st), K_fm=int(r.taper_size),
                      n_days=len(pmd.days_used))))
