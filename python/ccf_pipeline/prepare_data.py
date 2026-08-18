"""Translation of lib/ccf_prepare_data_Z.m, using obspy in place of the
hand-rolled load_sac.m/readsac.m SAC parsing (per direction: obspy handles
SAC I/O; what's translated here is the day-file pairing, validation,
distance filtering, and sliding-window cut-and-resample logic that obspy
doesn't provide for free).

Status: implemented but NOT yet verified against real data or against
Octave (no real SAC files were available while writing this -- see
NOTES.md). The windowing math (build_windows) is the part most worth
testing first since it's pure array logic, independent of file I/O.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from obspy import read, Trace
from obspy.geodetics import gps2dist_azimuth


@dataclass
class StationPairInfo:
    stanames: tuple[str, str]
    lats: tuple[float, float]
    lons: tuple[float, float]
    dt: float
    dist_km: float


def station_distance_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Interstation great-circle distance in km.

    CAVEAT: MATLAB's `distance()` (Mapping Toolbox) without an explicit
    ellipsoid argument computes angular distance on a *unit sphere*, and
    `deg2km()` converts using a fixed mean Earth radius (6371 km) -- a
    spherical-Earth approximation. obspy's `gps2dist_azimuth` uses the
    WGS84 *ellipsoid* by default (the modern standard). These differ by
    up to ~0.3% depending on latitude -- immaterial for a `dist_min`
    threshold filter at the km scale used here (20-40 km), but noted for
    fidelity. Not changed to the spherical formula since the ellipsoidal
    one is more accurate and the difference doesn't affect filtering
    outcomes at these thresholds.
    """
    dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return dist_m / 1000.0


def validate_pair(tr1: Trace, tr2: Trace, dt: float, min_samples: int = 30000) -> str | None:
    """Replicates ccf_prepare_data_Z.m's sanity checks. Returns None if the
    pair passes, or a short reason string if it should be skipped/errored.
    Matches the .m source's checks in the same order.
    """
    if tr1.stats.delta != tr2.stats.delta:
        raise ValueError("S1 and S2 sample rates don't match!")

    if abs(tr1.stats.starttime - tr2.stats.starttime) > tr1.stats.delta:
        raise ValueError("Station files do not have same start time")

    if abs(tr1.stats.delta - dt) >= 0.01 * dt or abs(tr2.stats.delta - dt) >= 0.01 * dt:
        raise ValueError("sampling interval does not match data! check dt")

    if np.all(tr1.data == 0) or np.all(tr2.data == 0):
        return "all zeros"

    if tr1.stats.npts == 0 or tr2.stats.npts == 0:
        return "no data"

    if tr1.stats.npts < min_samples:
        return f"sta1 too short ({tr1.stats.npts} samples)"
    if tr2.stats.npts < min_samples:
        return f"sta2 too short ({tr2.stats.npts} samples)"

    return None


def build_windows(
    tr1: Trace,
    tr2: Trace,
    winlength_hours: float,
    nstart_sec: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Core sliding-window cut-and-resample logic from ccf_prepare_data_Z.m
    (lines ~158-194), independent of file I/O. Returns (s1_windows,
    s2_windows), each shape (n_windows, win_length_samples) for this one
    day-pair -- the caller stacks these across days into the (day, window,
    sample) 3-D array the rest of the pipeline expects.

    Note: pts_begin/pts_end are computed exactly as in the .m source
    (1-based sample-index arithmetic, translated with a -1 offset at the
    point of array/interpolation-grid construction, not by changing the
    formula itself) to keep the window boundaries identical.
    """
    n1 = tr1.stats.npts
    nstart = nstart_sec / dt

    if winlength_hours == 24:
        nwin = 1
    else:
        nwin = int(np.floor(24 / winlength_hours)) * 2 - 1
        win_length = int(winlength_hours * 3600 / dt)
        last_pt = win_length * 0.5 * (nwin - 1) + 1 + nstart + win_length
        if last_pt < n1:
            nwin += 1

    win_length = int(winlength_hours * 3600 / dt)

    t1 = tr1.times()  # seconds from tr1.stats.starttime, i.e. S1Zt already zero-based
    t2_offset = tr2.stats.starttime - tr1.stats.starttime
    t2 = tr2.times() + t2_offset  # matches S2Zt = S2Zt + seconds(S2Ztstart-starttime)

    s1_windows = np.zeros((nwin, win_length + 1))
    s2_windows = np.zeros((nwin, win_length + 1))

    n2 = tr2.stats.npts
    for iwin in range(1, nwin + 1):  # 1-based to match the .m source's indexing formulas
        if winlength_hours == 24:
            pts_begin = nstart
            pts_end = n1 - nstart
        else:
            pts_begin = win_length * 0.5 * (iwin - 1) + 1 + nstart
            pts_end = pts_begin + win_length

        if pts_begin > n1 or pts_begin > n2 or pts_end > n1 or pts_end > n2:
            pts_begin = n1 - win_length - nstart
            pts_end = pts_begin + win_length

        tcut = np.arange(pts_begin, pts_end + 1) * dt  # +1: MATLAB's pts_begin:pts_end is inclusive

        s1_windows[iwin - 1, :] = np.nan_to_num(np.interp(tcut, t1, tr1.data, left=np.nan, right=np.nan))
        s2_windows[iwin - 1, :] = np.nan_to_num(np.interp(tcut, t2, tr2.data, left=np.nan, right=np.nan))

    return s1_windows, s2_windows


def find_day_pairs(datadir: Path, sta1: str, sta2: str, comp: str) -> list[tuple[Path, Path, str]]:
    """Lists sta1's day files and matches the corresponding sta2 file by
    filename suffix, replicating ccf_prepare_data_Z.m's file-discovery
    logic (lines ~48-70). Returns (path1, path2, day_id) triples for days
    where both stations have data.
    """
    dir1 = Path(datadir) / sta1
    pattern = f"*{comp}.sac"
    pairs = []
    for f1 in sorted(dir1.glob(pattern)):
        parts = f1.name.split(".")
        if len(parts) < 6:
            continue
        day_id = ".".join(parts[1:6])
        # sta2's file has the same suffix (everything after the station name)
        suffix = f1.name[len(sta1) :]
        f2_candidates = list((Path(datadir) / sta2).glob(f"{sta2}{suffix}"))
        if not f2_candidates:
            continue
        pairs.append((f1, f2_candidates[0], day_id))
    return pairs


def prepare_station_pair(
    datadir: Path,
    sta1: str,
    sta2: str,
    comp: str,
    winlength_hours: float,
    nstart_sec: float,
    dt: float,
    dist_min_km: float,
) -> tuple[np.ndarray, np.ndarray, StationPairInfo] | None:
    """Top-level equivalent of one station-pair's work inside
    ccf_prepare_data_Z.m's double loop. Returns (S1_data_mat, S2_data_mat,
    stapairsinfo) with the 3-D arrays shaped (day, window, sample), or
    None if the pair is skipped (too close, no overlapping days, etc.).
    """
    pairs = find_day_pairs(datadir, sta1, sta2, comp)
    if not pairs:
        return None

    s1_days = []
    s2_days = []
    pair_info = None

    for path1, path2, day_id in pairs:
        tr1 = read(str(path1))[0]
        tr2 = read(str(path2))[0]

        reason = validate_pair(tr1, tr2, dt)
        if reason is not None:
            continue

        if pair_info is None:
            lat1, lon1 = tr1.stats.sac.stla, tr1.stats.sac.stlo
            lat2, lon2 = tr2.stats.sac.stla, tr2.stats.sac.stlo
            dist_km = station_distance_km(lat1, lon1, lat2, lon2)
            if dist_km < dist_min_km:
                return None
            pair_info = StationPairInfo(
                stanames=(sta1, sta2), lats=(lat1, lat2), lons=(lon1, lon2), dt=dt, dist_km=dist_km
            )

        s1_w, s2_w = build_windows(tr1, tr2, winlength_hours, nstart_sec, dt)
        s1_days.append(s1_w)
        s2_days.append(s2_w)

    if not s1_days:
        return None

    return np.stack(s1_days, axis=0), np.stack(s2_days, axis=0), pair_info
