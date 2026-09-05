# `ADAMA_gvib.h5` verification: AF.SKRH-XV.BAND

Functional check of `ADAMA_gvib.h5` as a data source for the `findLowBand_ADAMAbenchmark`
notebook (`docs/findLowBand_ADAMAbenchmark_progress.md`): recover a known pair's raw traces from
`gvib.h5`, reproduce `ccf_prepare_data_T_mdg.m`'s windowing+rotation logic in Python, compute the
FastMspec cross-spectrum, and compare against the already-validated SKRH-BAND result from Stage 3
(`verification/skrh_band_real_data/`) -- checked by *reprocessing and comparing*, not by auditing
the SAC tree for completeness (that approach was tried, found unproductive, and abandoned per
direct guidance -- see the 2026-09-05 log in `docs/findLowBand_ADAMAbenchmark_progress.md`).

## Files

- `gvib_skrh_band_test.py`: loads `AF.SKRH`/`XV.BAND` raw Z/N/E traces directly from
  `ADAMA_gvib.h5` (via `h5py`, `obspyh5`'s own indexing convention), windows and rotates them to
  the Transverse component **entirely in memory** (no intermediate `.mat` file), and computes the
  FastMspec coherence. Saves the result to `.npz`.
- `compare_gvib_vs_reference.py`: three-way comparison -- (a) Sayan's own `matched_data.mat` run
  through our already-validated FastMspec, (b) the original precomputed MATLAB `coh_sum`
  (Stage 3's own reference), (c) this test's `gvib.h5`-sourced result.

## Simplifications relative to the full MATLAB script (stated explicitly, not hidden)

- **No instrument response removal**: confirmed both ways -- the actual production config
  (`a2_ccf_run_crosscorr_T_mdg.m`) has `IsRemoveIR=0`, and the SAC files themselves are already
  response-corrected (direct user confirmation) -- so using them as-is matches production exactly,
  not a shortcut.
- **Orientation correction assumed zero**: `OBS_orientations.txt` (a per-station H1-misalignment
  correction, needed for OBS deployments with unknown post-deployment orientation) was not
  located -- a tree-crawl search for it was abandoned per direct instruction rather than pursued
  further. `AF.SKRH` and `XV.BAND` are both standard land broadband stations, not OBS, so the
  correction this file provides very likely doesn't apply here -- a reasoned assumption, not a
  verified one.
- **Only full, clean 00:00:00-UTC-started day chunks are used**, skipping the MATLAB script's
  `interp1`-based handling of partial/odd-start-time recording boundaries. Reasonable for a
  validation test; not a production-grade port of the boundary logic.
- Channel naming (`.BHE`/`.BHN`/`.BHZ`, not `.BHR`/`.BHT`) confirms these are genuinely un-rotated
  components in `gvib.h5` -- the rotation step here is necessary, not redundant.

## A real bug found and fixed by this test (2026-09-05)

The first run produced a coherence that was the *same shape, opposite sign* of the known-good
reference at essentially every frequency -- caught by direct user inspection of the comparison
plot ("results are out of phase, suggesting errors in rotation..."), not by a passive metric.
Traced to the actual `rotate_vector.m`/`ccf_prepare_data_T_mdg.m` source, not re-derived from
memory: `ccf_prepare_data_T_mdg.m` rotates station 1 using `S1az` directly but station 2 using
`S2az + 180` -- the `+180` was missing from this script's station-2 rotation call. Separately
verified `rotate_vector.m`'s own trig formula (`vec_y = -sin(theta)*vec_1 + cos(theta)*vec_2`,
with the caller keeping `vec_y` as the transverse component) exactly matches this script's
`rotate_to_transverse`, so that part was already correct -- the missing `+180` was the only bug.
Fixed, rerun, and reconfirmed against the reference (see Result below). A sign-flip like this is
worth naming as a specific failure mode for whoever builds the eventual library: get one
station's azimuth convention backwards and the result still *looks* like real, structured
coherence (same shape) rather than obviously broken, making it easy to miss without a known-good
comparison to catch it against.

## Result (2026-09-05, after the rotation fix)

- **Distance cross-check**: 289.2 km (from `ADAMA_stalist.csv` coordinates, independent of Stage
  3's own SAC-header-derived value) vs. Stage 3's known 290.3 km -- 0.4% apart, confirms the right
  stations and geometry.
- **249 full-day, all-3-component overlapping days found** (2012-10-26 to 2013-07-01), all 15
  windows/day usable (`coh_num=3735 = 249*15` exactly, no truncation needed).
- **Sayan-sourced (our FastMspec on `matched_data.mat`) vs. original MATLAB `coh_sum`**: 2.7%
  relative L2 error -- reproduces Stage 3's own ~3% baseline, confirming this comparison's
  reference point is sound.
- **`gvib.h5`-sourced vs. both references, correctly in phase**: every peak and trough aligns with
  the reference across the entire 0-0.35 Hz band (the main burst at 0.06-0.11 Hz, secondary bursts
  at 0.16-0.20 and 0.25-0.31 Hz -- all match, not just approximately). Peak amplitude is roughly
  half the reference's (249 days here vs. 107 in Sayan's `matched_data.mat`) -- expected, not a
  bug: averaging over more independent days shrinks apparent coherence toward the true,
  lower-variance value, while a smaller sample shows inflated peaks. `matched_data.mat` is very
  likely an earlier, less complete snapshot of the same underlying archive.
- Plot: `docs/gvib_vs_reference_skrh_band.png`.

**Net conclusion**: `ADAMA_gvib.h5` is a trustworthy data source for this notebook, for at least
this pair, with the rotation now verified correct end-to-end against a known-good reference.
Doesn't rule out isolated failures elsewhere in the archive (per `ppToHDF5.py`'s silent-failure
risk, documented in the main tracker) -- this is one functional spot-check, not an exhaustive
guarantee.

## Performance/architecture note, for the eventual library (2026-09-05)

This test ran in a few minutes **on a login node**, not a dedicated compute job -- worth stating
why, since it bears directly on the library's design:
1. **No intermediate `.mat` round-trip.** Raw-trace read, windowing, rotation, and the FastMspec
   call all stay in memory as numpy arrays; no multi-GB array gets written to and read back from
   disk between steps.
2. **One open file handle, not thousands.** `gvib.h5`'s internal B-tree indexing lets `h5py` jump
   straight to each day/channel dataset; reading ~3,700 small datasets from one file avoids the
   per-file open/stat/close overhead a parallel filesystem pays for opening thousands of
   individual SAC files -- the architectural bet behind choosing `gvib.h5` over the SAC tree,
   confirmed empirically here, not just in principle.
3. **FastMspec itself is cheap by design** (`r/K` correction, no full-taper materialization) --
   the cross-spectrum computation over 3,735 windows was never the bottleneck.

**Memory footprint, measured, not estimated**: this one pair's full `S1_data_mat`/`S2_data_mat`
arrays are `(249, 15, 10801)` float64 each -- 323 MB per station, ~645 MB total for both. Small
enough that "load the whole pair into RAM at once" is the right default, not a premature
optimization -- but this won't hold for every pair (a pair with thousands of overlapping days, or
several pairs held concurrently by one worker, could exceed a shared node's per-task memory
budget). **Design principle for the eventual loading library, per direct user guidance**: memory-
aware by construction, not single-mode --
- Default: read a full pair's data into memory at once (as done here) when it comfortably fits.
- Fallback: chunk by day (or by however many days fit a given memory budget) and accumulate the
  coherence sum incrementally across chunks -- this is a natural fit for FastMspec's own
  cross-spectrum, which is fundamentally a sum/average across (day, window) pairs, so chunked
  accumulation doesn't require holding the full 3D array at once to get the same final result.
- Parallel-safe either way: `gvib.h5` supports concurrent read-only access from multiple `h5py`
  handles, so a SLURM array of workers, each assigned a different pair (or shard of pairs), can
  each open the file independently with no shared-state coordination needed.
Not yet implemented -- a design requirement to build the library against, not a retrofit.
