"""Pair-matched data loading from ADAMA_gvib.h5 -- the core asset for findLowBand_ADAMAbenchmark.

gvib.h5 is indexed per single-station-channel-day (obspyh5 convention:
`waveforms/{network}.{station}/{location}.{channel}/{starttime}_{endtime}`), not per-pair. This
module does the work of turning that into what compute_crosscorr_mtc_fastmspec (and the other
ccf_pipeline techniques) actually need: a matched, rotated, windowed (S1_data_mat, S2_data_mat)
pair, entirely from two station identifiers.

Validated against a real, already-known-good result (verification/gvib_skrh_band_test/,
AF.SKRH-XV.BAND): after fixing two real bugs found by that comparison --
(1) a rotation bug (station 2's azimuth needs +180 degrees, not just the raw back-azimuth --
confirmed directly against ccf_prepare_data_T_mdg.m/rotate_vector.m, not re-derived from memory),
and (2) a zero-filled-day dilution bug (AF.SKRH's data is 56.6% zero-filled on average across its
nominal "full" days -- a day-selection filter that only checks presence/length, not data validity,
silently lets these dead days dilute the coherence sum; fixed by porting
ccf_prepare_data_T_mdg.m's own unported "All zeros!" day-skip check, see
build_pair_matched_data's inline comment) -- this pipeline reproduces Stage 3's own validated
FastMspec/MATLAB result both in phase AND in amplitude across the full band (108 usable days after
exclusion, closely matching Sayan's own 107-day reference; coherence range -0.0920 to 0.0825 vs.
the reference's -0.1033 to 0.0834). The zero-day exclusion is a required correctness step for real
ADAMA stations, not optional hardening -- see verification/gvib_skrh_band_test/README.md for the
full hypothesis-elimination sequence (day-count/averaging, windowing-step-size, taper-bias,
normalization were each tested and ruled out before the true cause was found).

Design principles (per direct user guidance, 2026-09-05):
- **Memory-aware, not single-mode.** Measured: one pair's full arrays are ~650 MB for ~250
  overlapping days (323 MB per station, float64). Small enough that loading a full pair at once
  is the right default. `build_pair_matched_data`'s `chunk_days` parameter is the documented
  fallback for pairs/workloads that don't fit -- see its docstring for exactly what's implemented
  today and what's still a stub.
- **Parallel-safe by construction.** gvib.h5 supports concurrent read-only access from multiple
  `h5py.File` handles; nothing in this module holds cross-process state, so a SLURM array of
  workers can each open the file independently with no coordination needed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import h5py
import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth

DEFAULT_CHANNELS = ("BHZ", "BHN", "BHE")


def rotate_to_transverse(n: np.ndarray, e: np.ndarray, azimuth_deg: float) -> np.ndarray:
    """Standard N/E -> Transverse rotation, verified directly against rotate_vector.m
    (vec_y = -sin(theta)*vec_1 + cos(theta)*vec_2, caller keeps vec_y as the transverse
    component; vec_1=N, vec_2=E matches this project's own ch1/ch2 convention).

    Caller is responsible for the correct azimuth convention (see build_pair_matched_data's
    docstring) -- this function does not know which station it's rotating and cannot catch a
    missing +180 degrees itself. That mistake produces a coherence that is the *same shape,
    opposite sign* of the correct result -- easy to mistake for real structured coherence rather
    than obviously broken data (found and fixed exactly this way in
    verification/gvib_skrh_band_test/).
    """
    az = np.radians(azimuth_deg)
    return -n * np.sin(az) + e * np.cos(az)


def find_overlapping_full_days(f: h5py.File, sta1: str, sta2: str,
                                 channels=DEFAULT_CHANNELS) -> list[str]:
    """Dates (YYYY-MM-DD) where both stations have all requested channels present as a single
    chunk starting exactly at 00:00:00 UTC. Skips partial/odd-start-time recording-boundary
    chunks -- a real simplification (ccf_prepare_data_T_mdg.m instead interpolates across these
    boundaries), fine for a day-granularity pair match, not a full port of that boundary logic.
    """
    def day_set(sta, chan):
        grp = f["waveforms"][sta][f".{chan}"]
        return {k.split("_")[0].split("T")[0] for k in grp.keys()
                if k.split("_")[0].endswith("T00:00:00")}

    days = None
    for sta in (sta1, sta2):
        for chan in channels:
            d = day_set(sta, chan)
            days = d if days is None else (days & d)
    return sorted(days) if days else []


def load_day(f: h5py.File, sta: str, date: str, channels=DEFAULT_CHANNELS) -> dict[str, np.ndarray]:
    """Raw component arrays for one station on one full (00:00:00-started) day."""
    out = {}
    for chan in channels:
        grp = f["waveforms"][sta][f".{chan}"]
        key = next(k for k in grp.keys() if k.startswith(f"{date}T00:00:00"))
        out[chan] = grp[key][:]
    return out


@dataclass
class PairMatchedData:
    S1_data_mat: np.ndarray  # (n_day, n_window, n_samples)
    S2_data_mat: np.ndarray
    dist_km: float
    az_1to2_deg: float
    az_2to1_deg: float
    days_used: list = field(default_factory=list)


def build_pair_matched_data(gvib_path, sta1: str, sta2: str, stalist_path,
                              win_hours: float = 3.0, dt: float = 1.0, nstart_sec: float = 50.0,
                              channels=DEFAULT_CHANNELS,
                              chunk_days: int | None = None) -> PairMatchedData:
    """Main entry point: two station identifiers ("NET.STA") in, a ready-to-use
    (S1_data_mat, S2_data_mat) pair out -- matches ccf_prepare_data_T_mdg.m's own 50%-overlap
    windowing convention and rotation (verified, not re-derived -- see module docstring).

    Windowing, matched to ccf_prepare_data_T_mdg.m exactly, corrected here after being caught
    slightly wrong in the first version of this pipeline (verification/gvib_skrh_band_test/):
    the MATLAB "core" window length is `win_core = win_hours*3600/dt` samples (10800 for 3h,
    dt=1) -- this is what the 50%-overlap STEP uses (`win_core*0.5` per window), not the
    `win_core+1` sample COUNT each window actually contains (the "+1" comes only from MATLAB's
    inclusive-endpoint `pts_begin:pts_end` slicing). Using `win_core+1` for the step too (the
    original bug) drifts the window start by 0.5 samples per window, accumulating across a day's
    15 windows -- small in absolute terms but a real, avoidable deviation from exact fidelity.
    `nstart_sec` (`parameters.Nstart_sec` in `ccf_setup_params_T_mdg.m`, production value 50) is a
    per-day start offset in seconds, also omitted from the first version (assumed 0) -- included
    here as a parameter, defaulting to the confirmed production value.

    Rotation convention, stated explicitly since getting it wrong is easy to miss (see
    rotate_to_transverse's docstring): station 1's transverse uses the azimuth *from* sta1 *to*
    sta2 directly; station 2's uses the azimuth *from* sta2 *to* sta1 **plus 180 degrees**. This
    matches ccf_prepare_data_T_mdg.m exactly (its `S2az+180` term) and is confirmed correct by
    the AF.SKRH-XV.BAND validation (verification/gvib_skrh_band_test/): omitting the +180 gives a
    same-shape, opposite-sign result across the entire band.

    Orientation-file correction (OBS_orientations.txt, for non-standard station-frame offsets)
    is NOT applied here -- assumed zero, appropriate for standard land broadband stations, not
    yet confirmed for arbitrary ADAMA stations that might include OBS deployments. Check before
    trusting this function for a pair you haven't already validated by comparison, the way
    AF.SKRH-XV.BAND was.

    chunk_days: not yet implemented. When None (today's only working mode), loads every
    overlapping day into memory at once and returns the full 3D arrays -- fine for the ~650 MB /
    250-day pair this was validated against, but will not scale to a pair with a much larger day
    count or to a worker holding several pairs at once. The intended fallback (day-chunked
    incremental accumulation of the FastMspec coherence sum, since that sum is inherently
    separable across (day, window) pairs -- chunking changes memory profile, not the final
    numeric result) is a documented design target, not working code. Passing a non-None value
    raises NotImplementedError rather than silently doing the wrong thing.
    """
    if chunk_days is not None:
        raise NotImplementedError(
            "chunk_days is a documented design target, not yet implemented -- see this "
            "function's docstring. Passing None (the default) loads the full pair into memory."
        )

    stalist = pd.read_csv(stalist_path)

    def coords(net_sta):
        net, sta = net_sta.split(".")
        row = stalist[(stalist.Network == net) & (stalist.Station == sta)].iloc[0]
        return float(row.Latitude), float(row.Longitude)

    lat1, lon1 = coords(sta1)
    lat2, lon2 = coords(sta2)
    dist_m, az12, az21 = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    dist_km = dist_m / 1000.0

    win_core = int(win_hours * 3600 / dt)  # 10800 for 3h, dt=1 -- the MATLAB stepping unit
    win_length = win_core + 1  # 10801 -- the actual per-window SAMPLE COUNT (inclusive endpoint)
    nwin = int(24 // win_hours) * 2 - 1
    nstart = int(nstart_sec / dt)

    with h5py.File(gvib_path, "r") as f:
        days = find_overlapping_full_days(f, sta1, sta2, channels)
        if not days:
            raise ValueError(f"No overlapping full days found for {sta1}/{sta2}")

        s1_windows, s2_windows, used_days = [], [], []

        for date in days:
            d1 = load_day(f, sta1, date, channels)
            d2 = load_day(f, sta2, date, channels)
            minlen = min(len(d1["BHN"]), len(d1["BHE"]), len(d2["BHN"]), len(d2["BHE"]))
            if minlen < win_length + nstart:
                continue
            # ccf_prepare_data_T_mdg.m's own check, ported directly (was missing from the first
            # version of this pipeline -- found via a real ~40% amplitude discrepancy against
            # Sayan's own matched_data.mat, traced to this exact gap): a day that is entirely
            # zero (dead channel, gap filled with zeros rather than omitted -- confirmed present
            # for 249 of AF.SKRH's own "full" days at a 56.6% average zero-fraction, median
            # per-day RMS of exactly zero) still counts toward coh_num if not excluded, diluting
            # the coherence sum even though it contributes no real signal.
            if np.all(d1["BHN"][:minlen] == 0) or np.all(d2["BHN"][:minlen] == 0):
                continue
            t1 = rotate_to_transverse(d1["BHN"][:minlen], d1["BHE"][:minlen], az12)
            t2 = rotate_to_transverse(d2["BHN"][:minlen], d2["BHE"][:minlen], az21 + 180)

            day1 = np.zeros((nwin, win_length))
            day2 = np.zeros((nwin, win_length))
            for iwin in range(nwin):
                pts_begin = int(win_core * 0.5 * iwin) + nstart
                pts_end = pts_begin + win_length
                if pts_end > minlen:
                    pts_begin, pts_end = minlen - win_length, minlen
                day1[iwin, :] = t1[pts_begin:pts_end]
                day2[iwin, :] = t2[pts_begin:pts_end]
            s1_windows.append(day1)
            s2_windows.append(day2)
            used_days.append(date)

        if not used_days:
            raise ValueError(f"All overlapping days for {sta1}/{sta2} were excluded "
                              f"(too short or all-zero) -- no usable data.")

    return PairMatchedData(S1_data_mat=np.stack(s1_windows), S2_data_mat=np.stack(s2_windows),
                            dist_km=dist_km, az_1to2_deg=az12, az_2to1_deg=az21,
                            days_used=used_days)
