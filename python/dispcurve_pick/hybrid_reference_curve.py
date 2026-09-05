"""A general-purpose, per-path phase-velocity reference-curve library.

Built for Madagascar (this project's own 380 pairs) but designed to be called for *any* station
pair: pass two (lat, lon) points and a wave type ('love' or 'rayleigh'), get back a smooth,
physically-grounded reference curve for that specific great-circle path, instead of relying on one
generic curve for every path in a dataset. Directly targets the reference-curve-mismatch mechanism
documented in docs/round2_hypothesis_evaluation.tex Section 5.1.

Design goal (per direct user request): usable as a library by someone with no Madagascar-specific
context, and extensible toward global coverage without restructuring. Two things make that work:

1. **Sources are pluggable and geography-aware.** Each `PhaseVelocitySource` declares its own
   coverage region and period range and knows how to path-average itself; `build_reference_curve`
   just asks each configured source, in priority order, "do you cover this period at this
   location?" and uses whichever answers first. Today there are two:
   - `AdamaMap`: Africa only (short + long period, 6-40s here), from ADAMA_Maps
     (github.com/URseismology/ACE_ADAMA) -- the highest-resolution source where it applies.
   - `Gdm52Map`: global, but only >=25s (this module uses 45-150s to sit above ADAMA's range) --
     the universal fallback, everywhere ADAMA doesn't cover, and the only source at long period.
   Outside Africa, short periods (<25s) simply have no source yet and are left out of the curve
   rather than faked -- see "Roadmap" below.

2. **Every sample is tagged with its source and a `low_confidence` flag**, not just averaged
   silently into one curve. `AdamaMap`'s own short-period Love-wave uncertainty is well known (the
   user built ADAMA and is explicit that this is *part of why FastMSPEC exists* -- a more
   principled bandwidth choice is one route to doing better than ADAMA's own AkiEstimate pipeline
   at exactly the periods it's least confident in). Validated directly against a real ADAMA-measured
   pair (`XV.BITY-XV.MAPH`, see data/reference/hybrid_curve_README.md): good agreement at T=12-40s
   (-3.2% to +2.7%), much worse at T=6-10s (+15% to +95%, worst over an ocean-crossing path). Periods
   below `CAUTION_PERIOD_S` are flagged `low_confidence=True` so a caller (e.g. the picker) can
   widen its own search there rather than trust a single curve blindly.

Path-averaging convention (standard, Ekstrom-style): slowness-averaged, not velocity-averaged --
1/c_avg = mean(1/c(s)) over N points sampled uniformly along the great-circle path. This matches
how phase velocity actually combines along a path (each segment contributes travel time
proportional to its local slowness, not its local velocity).

Smoothing: raw per-period path averages are noisy (real tomographic-grid cell-to-cell uncertainty,
worse for short paths that average over few cells) -- fit with `scipy.interpolate.make_smoothing_spline`
(a proper penalized regression spline, GCV-selected smoothing), not an exact interpolant. Evaluated
against an unconstrained cubic `UnivariateSpline` first: that overshot to unphysical values (5.75
km/s against a raw data range of 3.8-4.8 km/s) near the noisy short-period end -- a real,
concretely-demonstrated risk, not a hypothetical one. `make_smoothing_spline` does not exhibit this
overshoot on the same data (see docs/smoothing_eval*.png).

Roadmap (per direct user direction, not yet built):
- Global short-period coverage: once FastMSPEC itself is validated at scale, its own large-scale
  measurements could become a third `PhaseVelocitySource` -- a FastMSPEC-derived global short-period
  map, plugged in ahead of GDM52 the same way ADAMA already is for Africa. The source-priority
  architecture above exists specifically so that addition doesn't require restructuring this module.
- Rayleigh: both ADAMA_Maps and GDM52 publish Rayleigh products in the same layout as Love
  (`RayAvgMap`/`R{period}_0_GDM52.pix`); `wave='rayleigh'` is wired through already, not yet
  validated against real data the way Love was.

Data provenance for both sources (not committed -- large/binary): data/reference/hybrid_curve_README.md.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.interpolate import make_smoothing_spline
from scipy.io import loadmat

CAUTION_PERIOD_S = 12.0
# Below this, treat the curve as low-confidence: validated (XV.BITY-XV.MAPH) short-period
# disagreement with real ADAMA measurements ran +15% to +95% for T=6-10s, vs. -3.2% to +2.7% for
# T=12-40s. Not a hard cutoff -- a caller-visible flag, since the true boundary is fuzzier than one
# number and depends on path (an ocean-crossing path was worse than an island-internal one at the
# same periods).

ADAMA_LOVE_PERIODS = [6, 8, 10, 12, 15, 20, 30, 35, 40]
ADAMA_RAYLEIGH_PERIODS = [6, 8, 10, 12, 15, 20, 25, 30, 35, 40]
# Love excludes 5s and 25s: ADAMA_Maps/LoveAvgMap only has _G_ (group velocity) files for those two
# periods, no _P_ (phase velocity) -- confirmed via the GitHub repo tree, not assumed. Using group
# velocity as a phase substitute would be physically wrong. Rayleigh's own availability has not
# been checked with the same rigor -- inherit the same list minus the two Love gaps as a starting
# assumption, flagged here for whoever validates the Rayleigh path.

GDM52_PERIODS = [45, 50, 60, 75, 100, 125, 150]

ADAMA_GRID_LAT_RANGE = (-33.92, 40.08)
ADAMA_GRID_LON_RANGE = (-23.20, 57.55)


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


class PhaseVelocitySource(ABC):
    """One period source (a regional or global tomographic map). Subclasses implement coverage
    and per-point lookup; path-averaging is shared here since the physics (slowness averaging)
    doesn't vary by source."""

    name: str

    @abstractmethod
    def periods_available(self, wave: str) -> list[float]:
        ...

    @abstractmethod
    def covers(self, lat: float, lon: float) -> bool:
        """Whether this source has (trustworthy) coverage at this point. Used to skip a source
        entirely for a path with any point outside its footprint, rather than silently averaging
        in extrapolated/garbage values from outside the grid's real coverage."""
        ...

    @abstractmethod
    def value_at(self, wave: str, period: float, lat: float, lon: float) -> float:
        ...

    def path_average_velocity(self, wave: str, period: float, lat1, lon1, lat2, lon2, n=50):
        lats, lons = great_circle_points(lat1, lon1, lat2, lon2, n)
        if not all(self.covers(la, lo) for la, lo in zip(lats, lons)):
            raise CoverageError(f"{self.name}: path leaves coverage region at T={period}s")
        vels = np.array([self.value_at(wave, period, la, lo) for la, lo in zip(lats, lons)])
        if np.any(~np.isfinite(vels)) or np.any(vels <= 0):
            raise CoverageError(f"{self.name}: invalid map values along path at T={period}s")
        slowness = 1.0 / vels
        return 1.0 / slowness.mean()


class CoverageError(Exception):
    """Raised when a source can't answer for a given path/period -- caller (build_reference_curve)
    catches this and moves to the next source in priority order."""


class AdamaMap(PhaseVelocitySource):
    """ADAMA_Maps (github.com/URseismology/ACE_ADAMA), Africa only, 0.25x0.25 degree grid.
    Highest-resolution source where it applies; also the one with known short-period uncertainty
    for Love waves (CAUTION_PERIOD_S)."""

    name = "ADAMA"

    def __init__(self, maps_dir: Path):
        self.maps_dir = Path(maps_dir)
        self._grids: dict[tuple[str, float], np.ndarray] = {}
        self.lat_axis = None
        self.lon_axis = None

    def periods_available(self, wave: str) -> list[float]:
        return ADAMA_LOVE_PERIODS if wave == "love" else ADAMA_RAYLEIGH_PERIODS

    def covers(self, lat: float, lon: float) -> bool:
        return (ADAMA_GRID_LAT_RANGE[0] <= lat <= ADAMA_GRID_LAT_RANGE[1] and
                ADAMA_GRID_LON_RANGE[0] <= lon <= ADAMA_GRID_LON_RANGE[1])

    def _prefix(self, wave: str) -> str:
        return "L" if wave == "love" else "R"

    def _subdir(self, wave: str) -> str:
        return "LoveAvgMap" if wave == "love" else "RayAvgMap"

    def _load(self, wave: str, period: float):
        key = (wave, period)
        if key in self._grids:
            return self._grids[key]
        prefix = self._prefix(wave)
        path = self.maps_dir / self._subdir(wave) / f"{prefix}{int(period)}_P_maps.mat"
        m = loadmat(path)
        avgmap = m["avgmap"]
        if self.lat_axis is None:
            latgrid = np.loadtxt(self.maps_dir / "a_latgrid_2Dgrid.txt", delimiter=",")
            longrid = np.loadtxt(self.maps_dir / "a_longrid_2Dgrid.txt", delimiter=",")
            self.lat_axis = latgrid[:, 0]
            self.lon_axis = longrid[0, :]
        self._grids[key] = avgmap
        return avgmap

    def value_at(self, wave: str, period: float, lat: float, lon: float) -> float:
        avgmap = self._load(wave, period)
        i = np.interp(lat, self.lat_axis, np.arange(len(self.lat_axis)))
        j = np.interp(lon, self.lon_axis, np.arange(len(self.lon_axis)))
        i0, j0 = int(np.floor(i)), int(np.floor(j))
        i1, j1 = min(i0 + 1, avgmap.shape[0] - 1), min(j0 + 1, avgmap.shape[1] - 1)
        di, dj = i - i0, j - j0
        v00, v01 = avgmap[i0, j0], avgmap[i0, j1]
        v10, v11 = avgmap[i1, j0], avgmap[i1, j1]
        return (v00 * (1 - di) * (1 - dj) + v01 * (1 - di) * dj +
                v10 * di * (1 - dj) + v11 * di * dj)


class Gdm52Map(PhaseVelocitySource):
    """Ekstrom's GDM52 (ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52), global, 1x1 degree
    geocentric grid, percent deviation from a per-period reference velocity (#PVEL0 header).
    The universal fallback source -- covers everywhere, but only >=25s."""

    name = "GDM52"

    def __init__(self, pix_dir: Path):
        self.pix_dir = Path(pix_dir)
        self._grids: dict[tuple[str, float], tuple] = {}

    def periods_available(self, wave: str) -> list[float]:
        return GDM52_PERIODS

    def covers(self, lat: float, lon: float) -> bool:
        return True  # global

    def _prefix(self, wave: str) -> str:
        return "L" if wave == "love" else "R"

    def _load(self, wave: str, period: float):
        key = (wave, period)
        if key in self._grids:
            return self._grids[key]
        prefix = self._prefix(wave)
        path = self.pix_dir / f"{prefix}{int(period):03d}_0_GDM52.pix"
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
        self._grids[key] = (pvel0, rows[:, 0], rows[:, 1], rows[:, 2])
        return self._grids[key]

    def value_at(self, wave: str, period: float, lat: float, lon: float) -> float:
        pvel0, lats, lons, devs = self._load(wave, period)
        # Geocentric-latitude correction not applied: negligible (<0.2 deg even at 45N/S) relative
        # to the 1-degree grid spacing -- a known simplification, not an oversight.
        lon180 = ((lon + 180) % 360) - 180
        idx = np.argmin((lats - lat) ** 2 + (lons - lon180) ** 2)
        return pvel0 * (1 + devs[idx] / 100.0)


@dataclass
class CurveSample:
    period: float
    velocity: float
    source: str
    low_confidence: bool


@dataclass
class ReferenceCurve:
    samples: list[CurveSample]
    func: object  # callable(freq_hz) -> velocity_km_s, smoothing-spline-based
    f_lo: float
    f_hi: float

    @property
    def periods(self):
        return np.array([s.period for s in self.samples])

    @property
    def velocities(self):
        return np.array([s.velocity for s in self.samples])

    @property
    def low_confidence_periods(self):
        return [s.period for s in self.samples if s.low_confidence]


def build_reference_curve(lat1, lon1, lat2, lon2, sources: list[PhaseVelocitySource],
                            wave: str = "love", n: int = 50,
                            caution_period_s: float = CAUTION_PERIOD_S) -> ReferenceCurve:
    """Build a smooth, per-path reference curve from a priority-ordered list of sources. Earlier
    sources in the list win for any period they cover; later sources (typically a lower-resolution
    global fallback) fill in periods/regions the earlier ones can't reach.

    wave: 'love' or 'rayleigh'.
    sources: e.g. [AdamaMap(...), Gdm52Map(...)] -- ADAMA's Africa-only, higher-resolution answer
        preferred where available, GDM52 filling in elsewhere and at long period. Order matters.
    """
    covered_periods: set[float] = set()
    samples: list[CurveSample] = []

    for source in sources:
        for period in source.periods_available(wave):
            if period in covered_periods:
                continue
            try:
                v = source.path_average_velocity(wave, period, lat1, lon1, lat2, lon2, n=n)
            except CoverageError:
                continue
            samples.append(CurveSample(period=period, velocity=v, source=source.name,
                                         low_confidence=period < caution_period_s))
            covered_periods.add(period)

    if len(samples) < 4:
        raise ValueError(f"Only {len(samples)} period samples available for this path -- "
                          f"too few to fit a reference curve (need >=4).")

    samples.sort(key=lambda s: s.period)
    periods = np.array([s.period for s in samples])
    vels = np.array([s.velocity for s in samples])
    freqs = 1.0 / periods
    order = np.argsort(freqs)
    freqs_sorted, vels_sorted = freqs[order], vels[order]

    spline = make_smoothing_spline(freqs_sorted, vels_sorted)
    return ReferenceCurve(samples=samples, func=spline, f_lo=freqs_sorted.min(),
                           f_hi=freqs_sorted.max())
