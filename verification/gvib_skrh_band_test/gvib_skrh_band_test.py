"""In-memory verification: recover AF.SKRH / XV.BAND raw traces directly from ADAMA_gvib.h5,
reproduce ccf_prepare_data_T_mdg.m's exact windowing+rotation logic in Python (no intermediate
.mat file -- everything held in numpy arrays), and compute the FastMspec cross-spectrum.

Compared against Stage 3's already-validated SKRH-BAND result (verification/skrh_band_real_data/):
if this pipeline (sourcing from gvib.h5 instead of Sayan's own matched_data.mat) reproduces a
consistent, physically sensible coherence, that's the functional check the user asked for --
NOT an inventory audit of the SAC tree, a reprocessing-and-compare check instead.

Simplifications relative to the full MATLAB script, stated explicitly:
- IsRemoveIR=0 in the actual production config (confirmed by reading a2_ccf_run_crosscorr_T_mdg.m
  directly) -- no instrument response removal needed; SAC data is already response-corrected
  (per direct user confirmation), used as-is either way.
- Orientation correction (OBS_orientations.txt) assumed zero: SKRH (AF) and BAND (XV) are both
  standard land broadband stations, not OBS -- the correction this file provides is specifically
  for OBS post-deployment orientation uncertainty, not applicable here. Not directly confirmed by
  reading that file (a tree-crawl search for it was abandoned per direct instruction); flagged as
  an assumption, not verified.
- Only full-day (00:00:00-started) chunks are used, skipping the MATLAB script's interpolation-
  based handling of partial/odd-start-time recording boundaries -- a reasonable simplification for
  a validation test, not a production-grade port.
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import sys
sys.path.insert(0, "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/python")

import numpy as np
import h5py
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
from datetime import datetime

from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec

GVIB_PATH = "/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5"
# Fetch: curl -sSL https://raw.githubusercontent.com/URseismology/ADAMA/main/DataFiles/ADAMA_stalist.csv
STALIST = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/ADAMA_stalist.csv"
WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5
DT = 1.0
WIN_HOURS = 3
WIN_LENGTH = int(WIN_HOURS * 3600 / DT) + 1  # 10801, matching this project's N_SAMPLES convention
NWIN = int(24 // WIN_HOURS) * 2 - 1  # 15, 50% overlap

STA1, STA2 = "AF.SKRH", "XV.BAND"


def rotate_to_transverse(n, e, azimuth_deg):
    """Standard 2D rotation: radial/transverse from N/E components, given the azimuth (degrees,
    clockwise from north) toward the other station. Matches rotate_vector's transverse output --
    T = -N*sin(az) + E*cos(az) is the standard back-azimuth-free transverse convention; verified
    against the sign convention implied by ccf_prepare_data_T_mdg.m's use for cross-correlation
    (only the relative sign between S1's and S2's transverse components matters for coherence
    magnitude; getting an overall sign flip would only affect the real-part sign, not convergence).
    """
    az = np.radians(azimuth_deg)
    return -n * np.sin(az) + e * np.cos(az)


def load_full_days(f, sta, chan):
    """Returns {date_str: (starttime_str, data_array)} for every full 24h-ish chunk starting at
    00:00:00 UTC for this station/channel."""
    grp = f["waveforms"][sta][f".{chan}"]
    out = {}
    for key in grp.keys():
        start_str, end_str = key.split("_")
        if not start_str.endswith("T00:00:00"):
            continue
        date = start_str.split("T")[0]
        out[date] = grp[key][:]
    return out


print(f"Loading traces for {STA1} and {STA2} from {GVIB_PATH} ...")
f = h5py.File(GVIB_PATH, "r")

days1 = {ch: load_full_days(f, STA1, ch) for ch in ["BHZ", "BHN", "BHE"]}
days2 = {ch: load_full_days(f, STA2, ch) for ch in ["BHZ", "BHN", "BHE"]}

common_days = set(days1["BHZ"]) & set(days1["BHN"]) & set(days1["BHE"]) & \
              set(days2["BHZ"]) & set(days2["BHN"]) & set(days2["BHE"])
common_days = sorted(common_days)
print(f"Full-day, all-3-component overlap: {len(common_days)} days "
      f"(range {common_days[0]} to {common_days[-1]})" if common_days else "NO OVERLAP FOUND")

# Station coordinates, from ADAMA's own station list (not re-derived from SAC headers here)
stalist = pd.read_csv(STALIST)
def coords(net_sta):
    net, sta = net_sta.split(".")
    row = stalist[(stalist.Network == net) & (stalist.Station == sta)].iloc[0]
    return float(row.Latitude), float(row.Longitude)

lat1, lon1 = coords(STA1)
lat2, lon2 = coords(STA2)
dist_m, az12, az21 = gps2dist_azimuth(lat1, lon1, lat2, lon2)
dist_km = dist_m / 1000.0
print(f"{STA1} ({lat1},{lon1}) -- {STA2} ({lat2},{lon2}): dist={dist_km:.1f} km, "
      f"az(1->2)={az12:.1f} deg, az(2->1)={az21:.1f} deg "
      f"(Stage 3's known value: 290.3 km -- cross-check)")

n_day = len(common_days)
S1_data_mat = np.zeros((n_day, NWIN, WIN_LENGTH))
S2_data_mat = np.zeros((n_day, NWIN, WIN_LENGTH))

for idate, date in enumerate(common_days):
    n1 = days1["BHN"][date]
    e1 = days1["BHE"][date]
    n2 = days2["BHN"][date]
    e2 = days2["BHE"][date]
    minlen = min(len(n1), len(e1), len(n2), len(e2))
    if minlen < WIN_LENGTH:
        continue
    # ccf_prepare_data_T_mdg.m: S1's transverse uses S1az directly; S2's uses S2az+180 (the
    # +180 was missing here originally -- a real bug, caught via out-of-phase comparison against
    # the known-good reference, not assumed correct).
    t1 = rotate_to_transverse(n1[:minlen], e1[:minlen], az12)
    t2 = rotate_to_transverse(n2[:minlen], e2[:minlen], az21 + 180)

    for iwin in range(NWIN):
        pts_begin = int(WIN_LENGTH * 0.5 * iwin)  # iwin=0-indexed here (MATLAB iwin-1)
        pts_end = pts_begin + WIN_LENGTH
        if pts_end > minlen:
            pts_begin = minlen - WIN_LENGTH
            pts_end = minlen
        S1_data_mat[idate, iwin, :] = t1[pts_begin:pts_end]
        S2_data_mat[idate, iwin, :] = t2[pts_begin:pts_end]

print(f"Assembled S1_data_mat/S2_data_mat: shape {S1_data_mat.shape}, in memory, no .mat file written")

r = compute_crosscorr_mtc_fastmspec(S1_data_mat, S2_data_mat, wband=WBAND, cutoff=CUTOFF, epsilon=EPSILON)
faxis = np.fft.fftfreq(WIN_LENGTH, d=DT)
pos = faxis > 0
coh = r.coh_sum[pos].real / r.coh_num
print(f"coh_num={r.coh_num}, coherence real-part range: {coh.min():.4f} to {coh.max():.4f}")

np.savez("/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/report_figures/gvib_skrh_band_coh.npz",
         freq=faxis[pos], coh=coh, dist_km=dist_km, n_day=n_day)
print("Saved coherence for comparison against Stage 3's known-good result.")
