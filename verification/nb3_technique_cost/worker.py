"""One technique, one bandwidth, in a fresh process -- so peak memory is a clean,
isolated number (ru_maxrss is a whole-process high-water mark; running techniques
back-to-back in one process lets later ones hide behind the first one's mark).

Usage: worker.py <s1_npy> <s2_npy> <technique> <nw> <out_npz> [k_taps]
  technique in {single-taper, MspecBestK, FastMspec, Mspec}
  k_taps (optional): force this taper count -- Mspec only, for the fixed-NW K-sweep
Writes <out_npz> with: freq, coh (real part, positive freqs), coh_num, taper_size
Prints one JSON line to stdout: {technique, nw, runtime_s, peak_mem_mb, taper_size, coh_num}
"""
import json, sys, time, resource
import numpy as np

sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")
from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import (
    compute_crosscorr_mtc_fastmspec, compute_crosscorr_mtc_mspec,
    compute_crosscorr_mtc_mspecbestk,
)

s1_npy, s2_npy, technique, nw_str, out_npz = sys.argv[1:6]
k_taps = int(sys.argv[6]) if len(sys.argv) > 6 else None
nw = float(nw_str)
s1 = np.load(s1_npy)
s2 = np.load(s2_npy)
if s1.ndim == 2:
    s1, s2 = s1[None], s2[None]
n_samples = s1.shape[2]
wband = nw / n_samples
CUTOFF, EPSILON = 1 - 1e-5, 1e-5

t0 = time.time()
if technique == "single-taper":
    # detrend + 5% cosine taper + plain FFT coherency -- work_unit.py's single-taper branch
    s1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1))
    s2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2))
    fft_s1, fft_s2 = np.fft.fft(s1p, axis=2), np.fft.fft(s2p, axis=2)
    coh_trace = fft_s2 * np.conj(fft_s1) / np.abs(fft_s1) / np.abs(fft_s2)
    coh_trace = np.where(np.isnan(coh_trace), 0, coh_trace)
    coh_sum = coh_trace.sum(axis=(0, 1))
    coh_num = coh_trace.shape[0] * coh_trace.shape[1]
    taper_size = 1
elif technique == "FastMspec":
    r = compute_crosscorr_mtc_fastmspec(s1, s2, wband=wband, cutoff=CUTOFF, epsilon=EPSILON)
    coh_sum, coh_num, taper_size = r.coh_sum, r.coh_num, r.taper_size
elif technique == "MspecBestK":
    r = compute_crosscorr_mtc_mspecbestk(s1, s2, wband=wband, cutoff=CUTOFF, epsilon=EPSILON, dt=1.0)
    coh_sum, coh_num, taper_size = r.coh_sum, r.coh_num, r.taper_size
elif technique == "Mspec":
    if k_taps is not None:
        # fixed-NW K-sweep (sec1b): vary K directly at a given NW
        r = compute_crosscorr_mtc_mspec(s1, s2, nw=nw, k_taps=k_taps, dt=1.0)
    else:
        # production config -- matches work_unit.py's batch pipeline exactly (NW_MSPEC=100, K_MSPEC=80),
        # the deliberately-heavy classical baseline behind every manifest / round2 number
        r = compute_crosscorr_mtc_mspec(s1, s2, nw=100, k_taps=80, dt=1.0)
    coh_sum, coh_num, taper_size = r.coh_sum, r.coh_num, r.taper_size
else:
    raise SystemExit(f"unknown technique {technique}")
runtime_s = time.time() - t0

faxis = np.fft.fftfreq(n_samples, d=1.0)
pos = faxis > 0
np.savez(out_npz, freq=faxis[pos], coh=(coh_sum[pos].real / coh_num),
         coh_num=coh_num, taper_size=taper_size, nw=nw, technique=technique)

peak_mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)  # Linux: KB -> MB... see note
# On Linux ru_maxrss is in kilobytes:
peak_mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
print(json.dumps(dict(technique=technique, nw=nw, runtime_s=round(runtime_s, 1),
                      peak_mem_mb=round(peak_mem_mb, 1), taper_size=int(taper_size),
                      coh_num=int(coh_num))))
