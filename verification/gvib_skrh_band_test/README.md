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

## A second real bug found and fixed (2026-09-05): zero-filled "dead" days diluting the average

After the rotation fix, `gvib.h5`-sourced coherence was correctly *in phase* with both references
across the whole band, but its peak amplitude was roughly half the reference's -- e.g. main-burst
peak ~0.04 vs. the reference's ~0.08. The first-pass explanation (below, now known wrong) was that
this was just an averaging-over-more-days effect: 249 overlapping full days found in `gvib.h5`
vs. 107 in Sayan's own `matched_data.mat`, and more independent days shrinks apparent coherence
toward a lower-variance true value. That explanation did not survive testing.

The user pushed back with a sequence of specific alternative hypotheses, each tested directly
rather than accepted or waved off:
- **Day-count/averaging dilution**: tested by randomly subsampling the 249-day set down to 107,
  50, and 20 days and recomputing. Amplitude stayed ~0.04-0.048 regardless of day count --
  ruled out. Real averaging dilution would show a clear day-count trend; none was found.
- **Windowing/overlap mismatch**: found and fixed a real (if separate) bug -- the 50%-overlap
  window *step* had used `win_length` (10801, the per-window sample count including MATLAB's
  inclusive-endpoint `+1`) instead of `win_core` (10800, the actual MATLAB stepping unit), and the
  per-day `Nstart_sec=50` production offset (`ccf_setup_params_T_mdg.m`) had been omitted entirely
  (assumed 0). Both fixed in `gvib_loader.py`. Amplitude was essentially unchanged by this fix
  (0.0398 -> 0.0399) -- a real correctness improvement, but not the cause of the discrepancy.
- **Taper bias / normalization differences**: reasoned through and ruled out -- taper bias would
  affect both day-counts identically (same `NW`/`K`), and FFT normalization was already controlled
  for since the Sayan-sourced run (same FastMspec code path) already matched the original MATLAB
  reference to ~3%.
- **Raw amplitude scale, checked directly per station**: compared RMS of `gvib.h5`'s raw rotated
  traces against Sayan's own `matched_data.mat` traces for the *same* days. Found the discrepancy
  was **station-specific, not global**: SKRH's `gvib.h5` amplitude was only ~0.60x Sayan's, while
  BAND's matched almost exactly (~1.0001x). A global effect (rotation, normalization, taper) would
  hit both stations equally -- this pointed straight at something wrong with SKRH's data
  specifically.
- **Zero-filled/dead-day gaps** (the user's own next hypothesis: "rotation angle orientation and
  possibly the new patch on Python code with zeros and nulls"): checked directly. **`AF.SKRH`'s
  `gvib.h5` data is 56.6% zero-filled on average across its nominal "full" days, with 141 of 249
  days having >1% zeros and a median per-day RMS of exactly zero** -- i.e. many of SKRH's "full,
  present" days are actually dead/gap-filled with zeros, not real recorded signal. `XV.BAND`, by
  contrast, was completely clean (0.0% zero-fraction on every day). The day-selection filter in
  the first version of `gvib_loader.py` checked only that a day chunk was *present* and *long
  enough* -- it never checked that the data inside was non-trivial. This is exactly the check
  `ccf_prepare_data_T_mdg.m` already has and skips explicitly (an `"All zeros!"` day-skip guard)
  that had not been ported into the Python loader.

**Root cause, confirmed**: roughly half of SKRH's "full days" in `gvib.h5` are zero-filled dead
data, silently diluting the coherence sum with (day, window) units that contribute zero real
signal but still count toward `coh_num` -- not a day-count/averaging effect, not a windowing bug,
not taper bias or normalization, and not a further rotation error (the rotation fix from the first
bug was already correct and stayed correct).

**Fix**: `gvib_loader.build_pair_matched_data` now excludes any day where either station's raw
N-component is entirely zero before windowing/rotation (matching MATLAB's own unported check),
switching the array assembly from a pre-sized fixed-day-count array to a dynamic list + `np.stack`
since the usable day count is now data-dependent.

## Result (2026-09-05, after both fixes)

- **Distance cross-check**: 289.2 km (from `ADAMA_stalist.csv` coordinates, independent of Stage
  3's own SAC-header-derived value) vs. Stage 3's known 290.3 km -- 0.4% apart, confirms the right
  stations and geometry.
- **108 usable days after zero-day exclusion** (down from 249 nominally-full days; `coh_num=1620
  = 108*15`), closely matching Sayan's own 107-day `matched_data.mat` -- strong independent
  evidence the fix identifies genuinely dead data rather than over-excluding.
- **Sayan-sourced (our FastMspec on `matched_data.mat`) vs. original MATLAB `coh_sum`**: 2.7%
  relative L2 error -- reproduces Stage 3's own ~3% baseline, confirming this comparison's
  reference point is sound.
- **Raw amplitude scale, corrected**: SKRH ratio improved from 0.60x to 0.911x Sayan's own values;
  BAND stayed at 0.994x (was 1.0001x) -- both now close to 1, consistent with the same underlying
  station data once dead days are excluded from both sides' effective sample.
- **`gvib.h5`-sourced vs. both references, in phase AND in amplitude**: coherence real-part range
  -0.0920 to 0.0825, closely matching the Sayan-sourced reference's -0.1033 to 0.0834 -- every peak
  and trough aligns across the entire 0-0.35 Hz band (the main burst at 0.06-0.11 Hz, secondary
  bursts at 0.16-0.20 and 0.25-0.31 Hz), and now the *heights* match too, not just the shape.
- Plot: `docs/gvib_vs_reference_skrh_band.png` (regenerated after the fix).

**Net conclusion**: `ADAMA_gvib.h5` is a trustworthy data source for this notebook, for at least
this pair, with both the rotation and the zero-day-exclusion now verified correct end-to-end
against a known-good reference. Doesn't rule out isolated failures elsewhere in the archive (per
`ppToHDF5.py`'s silent-failure risk, documented in the main tracker, and independently, the
zero-fill pattern found here for SKRH specifically) -- this is one functional spot-check, not an
exhaustive guarantee. **Practical implication for the eventual library/notebook**: any station may
have a nontrivial zero-filled fraction in its nominal "full" days; the zero-day-exclusion check
this bug forced into `gvib_loader.py` is not an edge-case nicety, it is necessary for correct
results on real ADAMA stations, and should be treated as a required step, not optional hardening.

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
