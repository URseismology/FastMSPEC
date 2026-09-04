# Hybrid per-pair reference curve: data provenance

`python/dispcurve_pick/hybrid_reference_curve.py` builds a per-pair Love-wave phase-velocity
reference curve from two path-averaged tomographic sources, replacing the single generic
`SDISPL.ASC` curve currently used for all 380 pairs.

## Sources

**ADAMA_Maps** (periods 6, 8, 10, 12, 15, 20, 30, 35, 40 s):
`github.com/URseismology/ACE_ADAMA`, `ADAMA_Maps/LoveAvgMap/L{period}_P_maps.mat` +
`ADAMA_Maps/LoveAvgMap/a_latgrid_2Dgrid.txt` / `a_longrid_2Dgrid.txt` for the grid axes (the
`.mat` files alone don't carry lat/lon; `plotpara.mat`, which the repo's own plotting scripts load
axes from, isn't checked into the repo -- these plain-text axis files are the substitute, found
by inspecting `ADAMA_Models/LoveMaps_txt/` instead).

**Real data gap, confirmed via the GitHub repo tree, not assumed**: periods 5s and 25s have only
`_G_` (group velocity) `.mat` files in `LoveAvgMap`, no `_P_` (phase velocity). Using group
velocity as a phase-velocity substitute would be physically wrong, so these two periods are simply
excluded from `ADAMA_PERIODS` rather than faked.

**GDM52** (periods 45, 50, 60, 75, 100, 125, 150 s): Ekström's global model,
`ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52`, `L{period}_0_GDM52.pix.gz`. Format: percent
deviation from a per-period reference velocity (`#PVEL0` header line), 1x1 degree geocentric grid
-- confirmed via the site's own `explain_GDM52_maps.txt`, not assumed from the filename. Absolute
velocity = `PVEL0 * (1 + deviation/100)`.

Neither dataset is committed to this repo (binary/large). Re-fetch with:
```bash
# ADAMA_Maps (per period)
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Maps/LoveAvgMap/L{period}_P_maps.mat" -o L{period}_P_maps.mat
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Models/LoveMaps_txt/a_latgrid_2Dgrid.txt" -o a_latgrid_2Dgrid.txt
curl -sSL "https://raw.githubusercontent.com/URseismology/ACE_ADAMA/main/ADAMA_Models/LoveMaps_txt/a_longrid_2Dgrid.txt" -o a_longrid_2Dgrid.txt

# GDM52 (per period)
curl -sSL "https://www.ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52/L{period:03d}_0_GDM52.pix.gz" -o L{period:03d}_0_GDM52.pix.gz
gunzip L{period:03d}_0_GDM52.pix.gz
```

## Method

Great-circle path between the two stations (spherical slerp, not ellipsoidal -- adequate given the
map resolutions of 0.25-1 degree already dominate achievable precision), sampled at 50 points.
Path-averaged phase velocity is **slowness-averaged** (standard, Ekström-style):
`1/c_avg = mean(1/c(s))` along the path, not a simple velocity average -- this matches how phase
velocity actually combines along a path (each segment contributes travel time proportional to its
local slowness).

## Validation (2026-09-04)

Tested against `XV.BITY-XV.MAPH`, a real ADAMA-measured pair (its own `cf` AkiEstimate dispersion
curve is known independently -- see `docs/notebook5_revamp_progress.md`'s 2026-09-04 log). Result:

- **T=12-40s: good agreement**, -3.2% to +2.7% (hybrid vs. ADAMA's real measured curve).
- **T=6-10s: large disagreement**, hybrid reads +15% to +95% too fast. This test pair crosses the
  open Mozambique channel; the tomographic map is likely poorly constrained at short periods over
  water (no seafloor stations feed the inversion there). Visually (not just by the percentages,
  which are inflated by the real curve's small absolute values at short period) the two curves
  still have a qualitatively similar rising shape, just offset -- see
  `docs/hybrid_vs_real_bity_maph.png`.
- A second check on `XV.BITY-XV.MAGY` (one of this project's own real, island-internal report
  example pairs, 223 km, no ocean crossing) does NOT show the same order-of-magnitude short-period
  blowup (values stay in a physically plausible 3.8-4.7 km/s range) -- but the curve is
  noticeably non-smooth/non-monotonic period-to-period, most likely real cell-to-cell uncertainty
  in the tomographic grid showing through because a short path averages over few grid cells (less
  noise cancellation than a long path). Not yet addressed with smoothing -- an open item.

**Net assessment**: promising for our actual use case (island-internal paths), but not yet blindly
production-ready. Two open items before wiring into the picker: (1) whether/how to smooth the
per-pair curve given real map-cell noise, (2) whether short periods (6-10s) need extra caution or
exclusion given the ocean-crossing test case's large disagreement, even though our own paths don't
cross ocean.
