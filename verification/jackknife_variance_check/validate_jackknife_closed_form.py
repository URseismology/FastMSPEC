"""Validate Haley & Anitescu (2017)'s Eq. 24 closed-form E{jackknife variance} against the
literal Eq. 23 delete-one jackknife, computed from REAL per-taper cross-spectra of the SKRH-BAND
pair (via Mspec's classical, per-taper-materializing approach -- unlike FastMspec, it keeps
individual eigenspectra, so it can serve as ground truth here). Run locally
(docs/notebook5_revamp_progress.md's Stage 5 log, 2026-09-03 entry) before trusting Eq. 24 as
FastMspec's variance proxy -- and it found a real gap: Eq. 24 assumes each per-taper term is
chi^2(2)-distributed, true for an auto-spectrum but not for a cross-spectrum term (the product
of two correlated complex Gaussians), whose distribution near a coherence null is qualitatively
different. See docs/stage5_bandwidth_theory.tex Section 5.1's 2026-09-03 update for the full
writeup this script's output grounds.

Data: same AFSKRH_XVBAND_win_3_all_matched_data.mat used by
verification/skrh_band_real_data/ -- see that folder's README for the retrieval command (not
committed here, matches this repo's established practice for bulk data). Place it alongside this
script, or pass its path as argv[1].
"""
import sys
from pathlib import Path
from math import comb

import numpy as np
from scipy.io import loadmat
from scipy.special import polygamma

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
from thomson_multitaper import dpss

PICK_FREQMIN, PICK_FREQMAX = 0.01, 0.5  # matches dispcurve_pick_batch/work_unit.py


def eq24_expected_jk_variance(K):
    """Haley & Anitescu supplement Eq. 24: E{jackknife variance}, a closed form in K alone."""
    trigamma = lambda t: polygamma(1, t)
    term1 = 2 / (K - 2) ** 2
    term2 = 0.5 * (trigamma((K - 1) / 2) - trigamma((K - 2) / 2))
    prefactor = (K - 1) ** 3 / (K * comb(K - 1, 2))
    return prefactor * (term1 + term2)


def literal_jk_variance_per_freq(per_taper_summed, K):
    """Eq. 23: literal delete-one jackknife variance of the log, per frequency bin.
    per_taper_summed: (n_freq, K) -- the K per-taper cross-spectrum magnitudes (already summed
    over traces/coh_num), kept individually rather than averaged over K.
    """
    total = per_taper_summed.sum(axis=1, keepdims=True)
    delete_one = (total - per_taper_summed) / (K - 1)
    log_delete_one = np.log(delete_one)
    log_delete_mean = log_delete_one.mean(axis=1, keepdims=True)
    return (K - 1) ** 2 / (K * comb(K - 1, 2)) * ((log_delete_one - log_delete_mean) ** 2).sum(axis=1)


def main():
    mat_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "AFSKRH_XVBAND_win_3_all_matched_data.mat"
    raw = loadmat(mat_path, simplify_cells=True)
    s1, s2 = raw["S1_data_mat"], raw["S2_data_mat"]
    if s1.ndim == 2:
        s1, s2 = s1[None, :, :], s2[None, :, :]
    n_day, n_window, n_samples = s1.shape
    n = n_samples
    x, y = s1.reshape(-1, n).T, s2.reshape(-1, n).T
    n_traces = x.shape[1]
    print(f"data shape (n_day, n_window, n_samples): {s1.shape}, n_traces={n_traces}, N={n}")

    faxis = np.fft.rfftfreq(n, d=1.0)
    band = (faxis >= PICK_FREQMIN) & (faxis <= PICK_FREQMAX)
    n_onesided = n // 2 + 1

    for K in [5, 13, 33]:
        NW = (K + 1) / 2  # rough inverse of K ~ 2NW-1
        tapers, _ = dpss(n, NW, K)

        # Per-taper cross-spectrum, kept individually -- exactly what FastMspec's architecture
        # deliberately avoids materializing (Section 4 of stage5_bandwidth_theory.tex), needed
        # here as ground truth for the literal jackknife.
        eigspec_k = np.zeros((n_onesided, n_traces, K), dtype=np.complex128)
        for k in range(K):
            fx = np.fft.fft(tapers[:, k : k + 1] * x, n=n, axis=0)[:n_onesided]
            fy = np.fft.fft(tapers[:, k : k + 1] * y, n=n, axis=0)[:n_onesided]
            eigspec_k[:, :, k] = fx * np.conj(fy)
        per_taper_summed = np.abs(eigspec_k.sum(axis=1))  # sum over traces (coh_num), keep K

        jk_var_literal = literal_jk_variance_per_freq(per_taper_summed, K)
        jk_var_closed_form = eq24_expected_jk_variance(K)

        mean_literal = jk_var_literal[band].mean()
        median_literal = np.median(jk_var_literal[band])
        print(f"\nK={K}: closed-form E{{sigma_J^2}} = {jk_var_closed_form:.6f}")
        print(f"       literal jackknife, band mean   = {mean_literal:.6f}")
        print(f"       literal jackknife, band median = {median_literal:.6f}")
        print(f"       ratio (mean literal / closed-form) = {mean_literal / jk_var_closed_form:.3f}")

        # Near-null diagnostic: bin by per-taper magnitude (proxy for proximity to a coherence
        # null) and check whether the discrepancy concentrates there.
        if K == 33:
            mag_band = per_taper_summed.mean(axis=1)[band]
            jkv_band = jk_var_literal[band]
            order = np.argsort(mag_band)
            q = len(order) // 4
            print("\n  Near-null diagnostic (K=33), jackknife variance by per-taper-magnitude quartile:")
            for i, label in enumerate(["lowest-mag (near-null-like)", "2nd", "3rd", "highest-mag"]):
                idx = order[i * q : (i + 1) * q] if i < 3 else order[3 * q :]
                print(f"    {label:28s}: mean={jkv_band[idx].mean():.4f}  median={np.median(jkv_band[idx]):.4f}")


if __name__ == "__main__":
    main()
