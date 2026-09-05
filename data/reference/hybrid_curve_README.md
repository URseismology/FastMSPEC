# Per-path reference curve library: data provenance

`python/dispcurve_pick/hybrid_reference_curve.py` is a general-purpose library: pass two (lat,
lon) points and a wave type, get a smooth, physically-grounded phase-velocity curve for that
specific great-circle path -- not tied to Madagascar, and designed for that path-independence to
extend cleanly toward global coverage. It replaces the single generic `SDISPL.ASC` curve currently
used for all 380 of this project's own pairs. Full design rationale is in the module's own
docstring; this file covers data provenance and the validation record.

## Sources (pluggable, priority-ordered)

**`AdamaMap`** -- ADAMA_Maps (`github.com/URseismology/ACE_ADAMA`), Africa only, 0.25x0.25 degree
grid, `ADAMA_Maps/{Love,Ray}AvgMap/{L,R}{period}_P_maps.mat` + `a_latgrid_2Dgrid.txt` /
`a_longrid_2Dgrid.txt` for the grid axes (the `.mat` files alone don't carry lat/lon;
`plotpara.mat`, which the repo's own plotting scripts load axes from, isn't checked into the repo
-- these plain-text axis files, found under `ADAMA_Models/LoveMaps_txt/`, are the substitute).
Love periods: 6, 8, 10, 12, 15, 20, 30, 35, 40 s.

**Real data gap, confirmed via the GitHub repo tree, not assumed**: periods 5s and 25s have only
`_G_` (group velocity) `.mat` files for Love, no `_P_` (phase velocity). Using group velocity as a
phase-velocity substitute would be physically wrong, so these two periods are simply excluded
rather than faked. (Rayleigh's own availability hasn't been checked with the same rigor --
`ADAMA_RAYLEIGH_PERIODS` currently just assumes the same gaps as Love, flagged in the module for
whoever validates the Rayleigh path.)

**`Gdm52Map`** -- Ekström's GDM52, global, 1x1 degree geocentric grid,
`ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52`, `{L,R}{period}_0_GDM52.pix.gz`. Format: percent
deviation from a per-period reference velocity (`#PVEL0` header line) -- confirmed via the site's
own `explain_GDM52_maps.txt`, not assumed from the filename. Absolute velocity =
`PVEL0 * (1 + deviation/100)`. Periods used here: 45, 50, 60, 75, 100, 125, 150 s (sits above
ADAMA's range; GDM52 itself starts at 25s but that's left to ADAMA where ADAMA covers it).

Neither dataset is committed to this repo (binary/large). Re-fetch with:
```bash
# ADAMA_Maps (per period, per wave prefix L or R)
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Maps/LoveAvgMap/L{period}_P_maps.mat" -o L{period}_P_maps.mat
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Models/LoveMaps_txt/a_latgrid_2Dgrid.txt" -o a_latgrid_2Dgrid.txt
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Models/LoveMaps_txt/a_longrid_2Dgrid.txt" -o a_longrid_2Dgrid.txt

# GDM52 (per period, per wave prefix L or R)
curl -sSL "https://www.ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52/L{period:03d}_0_GDM52.pix.gz" -o L{period:03d}_0_GDM52.pix.gz
gunzip L{period:03d}_0_GDM52.pix.gz
```
Expected local layout: `<adama_maps_dir>/LoveAvgMap/L{period}_P_maps.mat` +
`<adama_maps_dir>/a_latgrid_2Dgrid.txt`/`a_longrid_2Dgrid.txt`; `<gdm52_dir>/L{period:03d}_0_GDM52.pix`.

## Method

Great-circle path (spherical slerp, not ellipsoidal -- adequate given the map resolutions of
0.25-1 degree already dominate achievable precision), sampled at 50 points. Path-averaged phase
velocity is **slowness-averaged** (standard, Ekström-style): `1/c_avg = mean(1/c(s))` along the
path, not a simple velocity average -- matches how phase velocity actually combines along a path.

**Smoothing**: raw per-period path averages are noisy (real tomographic-grid cell-to-cell
uncertainty, worse for short paths averaging over few cells). Evaluated an unconstrained cubic
`scipy.interpolate.UnivariateSpline` first -- it overshot to unphysical values (5.75 km/s against
a raw data range of 3.8-4.8 km/s) right at the noisy short-period end, a concretely demonstrated
risk (`docs/smoothing_eval.png`). Switched to `scipy.interpolate.make_smoothing_spline` (a proper
penalized regression spline, GCV-selected smoothing) -- stays within the data range, no overshoot,
smooth without chasing individual-cell noise (`docs/smoothing_eval2.png`).

## Validation (2026-09-04)

Tested against `XV.BITY-XV.MAPH`, a real ADAMA-measured pair (its own `cf` AkiEstimate dispersion
curve is known independently -- see `docs/notebook5_revamp_progress.md`'s 2026-09-04 log).

- **T=12-40s: good agreement**, -3.2% to +2.7% (hybrid vs. ADAMA's real measured curve).
- **T=6-10s: large disagreement**, hybrid reads +15% to +95% too fast. This test pair crosses the
  open Mozambique channel; the tomographic map is likely poorly constrained at short periods over
  water (no seafloor stations feed the inversion there). Visually the two curves still have a
  qualitatively similar rising shape, just offset -- the alarming percentages are partly an
  artifact of the real curve's small absolute values at short period (`docs/hybrid_vs_real_bity_maph.png`).
- A second check on `XV.BITY-XV.MAGY` (this project's own island-internal report example, 223 km,
  no ocean crossing) does not show the same magnitude of short-period blowup (stays physically
  plausible, 3.8-4.7 km/s) but is noticeably noisy period-to-period before smoothing
  (`docs/hybrid_bity_magy_island.png`) -- addressed by the smoothing-spline switch above.

**Context that matters for how to read this** (direct user note): ADAMA's own short-period Love
uncertainty is a *known, expected* limitation -- the user built ADAMA and states this uncertainty
is part of the motivation for FastMSPEC itself (a more principled bandwidth choice is one route to
doing better than AkiEstimate at exactly the periods it's least confident in). This isn't a defect
in the reference-curve approach to quietly work around; it's the reason this whole project exists,
and it's handled explicitly, not hidden: every `CurveSample` below `CAUTION_PERIOD_S` (12s, set
from this validation) carries `low_confidence=True` so a caller can react accordingly, rather than
trusting the whole curve uniformly.

## Status

Library-ready (pluggable `PhaseVelocitySource` sources, wave-type parameter, per-sample
source/confidence tagging) but **not yet wired into `python/dispcurve_pick_batch/work_unit.py`**
-- that integration, plus what "extra caution below 12s" should concretely mean for the picker
(e.g. a period-dependent corridor width in `build_template_family`, vs. down-weighting
low-confidence periods in scoring), is the next step, not yet decided or implemented.

## Roadmap (per direct user direction, not yet built)

- **Global short-period source**: once FastMSPEC is validated at scale, its own large-scale
  measurements could become a third `PhaseVelocitySource` -- a FastMSPEC-derived global
  short-period map, positioned ahead of GDM52 (and, outside Africa, filling the gap ADAMA leaves)
  in the same priority-ordered way `AdamaMap` already sits ahead of `Gdm52Map`. This is exactly
  why sources are pluggable rather than hardcoded to two: adding a third shouldn't require
  restructuring this module.
- **Rayleigh**: `wave='rayleigh'` is wired through the API already (both ADAMA_Maps and GDM52
  publish Rayleigh products in the same layout as Love) but not validated against real data the
  way Love was here.
