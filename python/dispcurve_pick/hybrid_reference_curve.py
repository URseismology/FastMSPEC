"""Hybrid per-pair Love-wave phase-velocity reference curve: ADAMA_Maps (5-40s) + GDM52 (>40s),
path-averaged along each station pair's great-circle path.

Replaces the single generic SDISPL.ASC curve family (one curve used for all 380 pairs) with a
physically-grounded, per-pair curve -- directly targeting the reference-curve-mismatch mechanism
documented in docs/round2_hypothesis_evaluation.tex Section 5.1.

Data provenance (not committed -- see README alongside this module):
- ADAMA_Maps: github.com/URseismology/ACE_ADAMA, ADAMA_Maps/LoveAvgMap/L{period}_P_maps.mat.
  0.25x0.25 degree grid, lat -33.92..40.08 (row 0 = south), lon -23.20..57.55 (col 0 = west).
  Periods 6,8,10,12,15,20,30,35,40 s (5s and 25s excluded -- see ADAMA_PERIODS below). Values in
  km/s, no ocean mask applied (raw .mat) -- the mask matters for African-continent-wide use but
  not for Madagascar-internal paths, which stay on land.
- GDM52: ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52, L{period}_0_GDM52.pix.gz. 1x1 degree
  geocentric grid, percent deviation from a per-period reference velocity (#PVEL0 header line).
  Periods used here: 45,50,60,75,100,125,150 s (the >40s extension).

Path-averaging convention (standard, Ekstrom-style): slowness-averaged, not velocity-averaged --
1/c_avg = mean(1/c(s)) along N points sampled uniformly on the great-circle path. This matches
how phase velocity actually combines along a path (a group of N points each contribute travel
time proportional to their local slowness).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat

ADAMA_PERIODS = [6, 8, 10, 12, 15, 20, 30, 35, 40]
# Periods 5s and 25s are excluded: ADAMA_Maps/LoveAvgMap only has _G_ (group velocity) files for
# these two periods, no _P_ (phase velocity) -- confirmed via the GitHub repo tree, not assumed.
# Using group velocity as a phase-velocity substitute would be physically wrong, so these two
# periods are simply missing from the reference curve rather than faked.
GDM52_PERIODS = [45, 50, 60, 75, 100, 125, 150]
EARTH_RADIUS_KM = 6371.0


def great_circle_points(lat1, lon1, lat2, lon2, n=50):
    """N points uniformly spaced (by fraction of angular distance) along the great-circle path,
    via spherical slerp. Spherical approximation, not ellipsoidal -- adequate given the map
    resolutions (0.25-1 degree) already dominate the achievable precision here.
    """
    p1 = np.radians([lat1, lon1])
    p2 = np.radians([lat2, lon2])
    x1 = np.array([np.cos(p1[0]) * np.cos(p1[1]), np.cos(p1[0]) * np.sin(p1[1]), np.sin(p1[0])])
    x2 = np.array([np.cos(p2[0]) * np.cos(p2[1]), np.cos(p2[0]) * np.sin(p2[1]), np.sin(p2[0])])
    dot = np.clip(np.dot(x1, x2), -1.0, 1.0)
    omega = np.arccos(dot)
    if omega < 1e-10:
        return np.full(n, lat1), np.full(n, lon1)
    t = np.linspace(0, 1, n)
    s1 = np.sin((1 - t) * omega) / np.sin(omega)
    s2 = np.sin(t * omega) / np.sin(omega)
    x = np.outer(s1, x1) + np.outer(s2, x2)
    lat = np.degrees(np.arcsin(np.clip(x[:, 2], -1, 1)))
    lon = np.degrees(np.arctan2(x[:, 1], x[:, 0]))
    return lat, lon


class AdamaLoveMap:
    """Reads ADAMA_Maps/LoveAvgMap/L{period}_P_maps.mat and provides bilinear lookup."""

    def __init__(self, maps_dir: Path):
        self.maps_dir = Path(maps_dir)
        self._grids = {}  # period -> (lat_axis, lon_axis, avgmap)
        self.lat_axis = None
        self.lon_axis = None

    def _load(self, period):
        if period in self._grids:
            return self._grids[period]
        path = self.maps_dir / f"L{period}_P_maps.mat"
        m = loadmat(path)
        avgmap = m["avgmap"]
        if self.lat_axis is None:
            lat_path = self.maps_dir / "a_latgrid_2Dgrid.txt"
            lon_path = self.maps_dir / "a_longrid_2Dgrid.txt"
            latgrid = np.loadtxt(lat_path, delimiter=",")
            longrid = np.loadtxt(lon_path, delimiter=",")
            self.lat_axis = latgrid[:, 0]
            self.lon_axis = longrid[0, :]
        self._grids[period] = avgmap
        return avgmap

    def value_at(self, period, lat, lon):
        avgmap = self._load(period)
        # bilinear interpolation on the regular grid
        i = np.interp(lat, self.lat_axis, np.arange(len(self.lat_axis)))
        j = np.interp(lon, self.lon_axis, np.arange(len(self.lon_axis)))
        i0, j0 = int(np.floor(i)), int(np.floor(j))
        i1, j1 = min(i0 + 1, avgmap.shape[0] - 1), min(j0 + 1, avgmap.shape[1] - 1)
        di, dj = i - i0, j - j0
        v00, v01 = avgmap[i0, j0], avgmap[i0, j1]
        v10, v11 = avgmap[i1, j0], avgmap[i1, j1]
        return (v00 * (1 - di) * (1 - dj) + v01 * (1 - di) * dj +
                v10 * di * (1 - dj) + v11 * di * dj)

    def path_average_velocity(self, period, lat1, lon1, lat2, lon2, n=50):
        lats, lons = great_circle_points(lat1, lon1, lat2, lon2, n)
        vels = np.array([self.value_at(period, la, lo) for la, lo in zip(lats, lons)])
        if np.any(~np.isfinite(vels)) or np.any(vels <= 0):
            raise ValueError(f"ADAMA map has invalid values along path at T={period}s "
                              f"({np.sum(~np.isfinite(vels) | (vels <= 0))} of {n} points)")
        slowness = 1.0 / vels
        return 1.0 / slowness.mean()


class Gdm52LoveMap:
    """Reads GDM52 L{period}_0_GDM52.pix files (percent deviation from PVEL0) and provides
    nearest-cell lookup on the 1-degree geocentric grid."""

    def __init__(self, pix_dir: Path):
        self.pix_dir = Path(pix_dir)
        self._grids = {}  # period -> (pvel0, lat_array, lon_array, dev_array)

    def _load(self, period):
        if period in self._grids:
            return self._grids[period]
        path = self.pix_dir / f"L{period:03d}_0_GDM52.pix"
        pvel0 = None
        rows = []
        with open(path) as f:
            for line in f:
                if line.startswith("#PVEL0"):
                    pvel0 = float(line.split(":")[1])
                elif not line.startswith("#"):
                    lat, lon, pix, dev = map(float, line.split())
                    rows.append((lat, lon, dev))
        rows = np.array(rows)
        self._grids[period] = (pvel0, rows[:, 0], rows[:, 1], rows[:, 2])
        return self._grids[period]

    def value_at(self, period, lat, lon):
        pvel0, lats, lons, devs = self._load(period)
        # geocentric latitude correction: negligible relative to 1-degree grid spacing at these
        # latitudes (<0.2 deg even at 45N/S) -- not applied, noted as a known simplification.
        lon180 = ((lon + 180) % 360) - 180
        idx = np.argmin((lats - lat) ** 2 + (lons - lon180) ** 2)
        return pvel0 * (1 + devs[idx] / 100.0)

    def path_average_velocity(self, period, lat1, lon1, lat2, lon2, n=50):
        lats, lons = great_circle_points(lat1, lon1, lat2, lon2, n)
        vels = np.array([self.value_at(period, la, lo) for la, lo in zip(lats, lons)])
        slowness = 1.0 / vels
        return 1.0 / slowness.mean()


@dataclass
class HybridCurveResult:
    periods: np.ndarray
    velocities: np.ndarray
    func: object  # interp1d, velocity(freq_hz)
    f_lo: float
    f_hi: float


def build_hybrid_reference_curve(lat1, lon1, lat2, lon2, adama_maps_dir, gdm52_dir,
                                   adama_periods=ADAMA_PERIODS, gdm52_periods=GDM52_PERIODS, n=50):
    """Path-averaged hybrid Love-wave reference curve for one station pair. Returns a
    HybridCurveResult with .func(freq_hz) -> velocity_km_s, matching load_reference_curve's
    (c_ref, f_lo, f_hi) interface when used as (result.func, 1/max(periods), 1/min(periods)).
    """
    adama = AdamaLoveMap(adama_maps_dir)
    gdm52 = Gdm52LoveMap(gdm52_dir)

    periods, vels = [], []
    for T in adama_periods:
        v = adama.path_average_velocity(T, lat1, lon1, lat2, lon2, n=n)
        periods.append(T)
        vels.append(v)
    for T in gdm52_periods:
        v = gdm52.path_average_velocity(T, lat1, lon1, lat2, lon2, n=n)
        periods.append(T)
        vels.append(v)

    periods = np.array(periods, dtype=float)
    vels = np.array(vels, dtype=float)
    freqs = 1.0 / periods  # descending as periods ascend
    order = np.argsort(freqs)
    freqs, vels_sorted = freqs[order], vels[order]

    func = interp1d(freqs, vels_sorted, bounds_error=False, fill_value=(vels_sorted[0], vels_sorted[-1]))
    return HybridCurveResult(periods=periods, velocities=vels, func=func,
                              f_lo=freqs.min(), f_hi=freqs.max())
