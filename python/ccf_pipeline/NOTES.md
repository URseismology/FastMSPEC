# Translation notes: `ccf_compute_crosscorr_mtc_Z.m` -> `ccf_pipeline`

[← Back to repo README](../../README.md) | See also: [verification/octave_verify_ccf_pipeline/README.md](../../verification/octave_verify_ccf_pipeline/README.md) | [legacy/matlab_source/README.md](../../legacy/matlab_source/README.md)

See [`docs/plan_ccf_mtc_translation.md`](../../docs/plan_ccf_mtc_translation.md) (repo root) for the full phased plan and the resolved
call-graph investigation. This file tracks per-phase implementation status and caveats, same
role as [`python/NOTES.md`](../NOTES.md) for the multitaper library.

## Status

**Phase 0 (done)**: empirically resolved the axis-semantics question (no bug -- `N` is
genuinely the sample count) and found a real bug in the one-sided-to-two-sided reflection step.
Full writeup: `../../verification/octave_verify_ccf_pipeline/README.md`.

**Phase 1 (done)**: `'FastMspec'` technique for the Z component
(`crosscorr_mtc.compute_crosscorr_mtc_fastmspec`). Verified against a real Octave run of the
unmodified `ccf_compute_crosscorr_mtc_Z.m` (`tests/test_crosscorr_mtc.py`, 4/4 passing):
- The one-sided spectrum (unaffected by the reflection bug) matches to relative L2 error ~9e-16.
- `coh_num` and `taper_size` (K + transition-region tapers) match exactly.
- The two-sided output length is correct (32, matching `N`) where Octave's literal buggy output
  is 33 -- the fixture (`N=32`, `N mod 4 == 0`) was deliberately chosen to land in the bug's
  trigger zone, so this test doubles as a regression check for the fix.

**Phase 2 (done)**: `'Mspec'` and `'MspecBestK'` techniques
(`compute_crosscorr_mtc_mspec`, `compute_crosscorr_mtc_mspecbestk`), using the existing `dpss()`
wrapper in place of `sleptap()` and a new `classical_spectrum_batch` (translating `mspec_fast.m`'s
plain `avgspec` subfunction -- simple taper-averaged cross-periodogram, no sinc-kernel fusion).
Both verified against real Octave runs on the same fixture as Phase 1: one-sided spectrum matches
to ~2e-15 relative error, taper counts match exactly (3 for `Mspec`, 2 for `MspecBestK`). 6/6 tests
passing.

**Found a second real upstream bug while verifying `MspecBestK`**: the original MATLAB branch
calls `mspec_fast` requesting only 4 output arguments (no `totalMB`), then unconditionally reads
`saved_ccf_path.psi_memory_space = totalMB` afterward -- this is an undefined-variable error at
runtime, confirmed by actually running the unmodified `.m` file in Octave (it crashes exactly
there). The coherency computation and its `save()` call both happen *before* this crash, so a
valid reference `coh_sum` was still recoverable for verification. The Python port implements the
(correct, verifiable) coherency computation and simply omits the memory-diagnostic field --
not a translation gap, the original crashes here every time this branch runs.

**Phase 3 (done)**: `_T.m` sibling's calling-convention differences, supported via optional
parameters on the same Z functions rather than duplicated code (`compute_crosscorr_mtc_fastmspec`
accepts an optional pre-built `fmtse=`, matching `_T.m`'s `filttype.FMTSE` injection pattern;
`compute_crosscorr_mtc_mspec` accepts optional direct `nw`/`k_taps`, matching `_T.m`'s
`NW_mspec`/`K_taps_mspec` config fields). `MspecBestK` is identical between Z and T, no change
needed. Both new calling conventions verified against real Octave runs of the actual
`ccf_compute_crosscorr_mtc_T.m` (not just `_Z.m`): ~1e-15 relative error. 8/8 tests passing at
this point.

**Phase 4 (mostly done)**: preprocessing stage (`preprocessing.py`) -- `ccf_detrend_3dim` and
`ccf_cos_taper_3dim` verified against real Octave runs to machine precision (~1e-12 and exact,
respectively). `ccf_butterfilt_3dim` runs and produces finite, correctly-shaped output, but its
exact match against the real `FiltFiltM.m` was **not achieved** -- see the "Known unresolved
discrepancy" section below. Deliberately scoped to only detrend/taper/prefilter: the
`IsMultiTaper`/`IsFTN`/`IsOBN`/`IsSpecWhiten` stages in `ccf_preprocess_filter_data.m` feed a
*different* downstream pathway (the outer dispatcher's already-frequency-domain branch) that the
real production config found in this codebase doesn't use alongside `IsMspec=1` -- see
`preprocessing.py`'s module docstring.

**Phase 5 (done, real-data verified -- see "Real-data verification" section below)**: obspy-based SAC loading/windowing
(`prepare_data.py`), translating `ccf_prepare_data_Z.m`'s day-file pairing, validation checks,
min-distance filter, and sliding-window cut-and-resample logic (obspy's `read()` handles the SAC
parsing itself, so `load_sac.m`/`readsac.m` weren't translated). No real SAC files were available
while writing this, so verification so far is self-consistency only:
- `build_windows`' output shape matches the `.m` source's own `nwin`/`win_length` formula exactly.
- Window boundaries land on the exact expected sample indices (checked with a ramp signal where
  window contents directly reveal which indices were cut).
- `validate_pair` correctly accepts good synthetic data and rejects all-zero/too-short traces,
  matching the `.m` source's checks and their order.

7/7 self-consistency tests passing; **since updated with real-data/Octave verification, see
"Real-data verification" below**. Distance calculation uses obspy's WGS84-ellipsoid
`gps2dist_azimuth` rather than MATLAB's spherical-Earth `distance()`+`deg2km()` -- documented in
`station_distance_km`'s docstring as an intentional, immaterial-at-these-thresholds difference,
not an oversight.

**Phase 6 (dispatcher done, real-data verified for the plain-fft branch -- see "Real-data verification" below)**: `dispatch.py`'s `compute_crosscorr`
translates `ccf_compute_crosscorr_Z.m`/`_T.m`'s three-way dispatcher (`IsMultiTaper||IsFTN||IsOBN
||IsSpecWhiten` / `IsMspec` / plain-fft). The `IsMspec` routing -- the actual production path --
is verified against a real Octave run of the **dispatcher itself** (not just `_mtc_Z.m` called
directly), confirming the wiring is correct end-to-end: ~9e-16 relative error, exact taper/trace
counts. The other two branches are implemented for fidelity to the source's structure and
smoke-tested (finite output, correct shape) but are not the exercised production path (every real
config found sets `IsMspec=1`) and have no Octave cross-check. Multi-level stacking
(`ccf_save_computed_ccf_Z`'s day/month/single-stack aggregation) is **not translated** -- every
real config found only uses `IsOutputFullstack=1`, which the `_mtc_` functions already handle via
their own internal full-stack summation.

**Update**: the real end-to-end run described above as a gap has since been done -- see
"Real-data verification" below. What that section does *not* cover: the `IsMspec`/`FastMspec`
path (Phases 0-3) against real data (only the plain-fft branch was real-data-verified so far),
and T-component rotation against real N/E data.

**Test suite status**: 23/23 tests passing across the whole `ccf_pipeline` package
(`tridieig`/`tridisolve`-level checks aside, which live in the separate `thomson_multitaper`
package).

## Known unresolved discrepancy: `ccf_butterfilt_3dim` / `FiltFiltM.m`

`FiltFiltM.m`'s own header says it implements Gustafsson's (1996) zero-phase initial-condition
method, "rewritten from scratch" (its own changelog). Tried both scipy `filtfilt` defaults
(reflect-padding) and `method='gust'` (scipy's own Gustafsson implementation) against a real
Octave run of the actual `ccf_butterfilt_3dim.m` + `FiltFiltM.m` (once `octave-signal` was
installed via `sudo apt install octave-signal`, needed for `butter()`):
- Short (32-sample) test window: up to 19-50% relative error depending on method, worst near a
  few specific traces.
- Longer (500-sample) test window, excluding a 20-sample margin from each edge: still ~1-2%
  relative error -- too large to be edge-transient noise, and not shrinking toward zero the way
  a pure edge/padding-convention mismatch would.

This looks like a genuine difference in how `FiltFiltM.m`'s specific (non-standard) Gustafsson
matrix formulation solves for initial state, not just a boundary-padding choice. Not pursued
further: `filttype.IsPrefilter = 0` in every real production config found in this codebase (both
`a2_ccf_run_crosscorr_T_mdg.m` and `ccf_synth_setup_perfrm_evl.m`), so this function is not on
the path that's actually exercised. If prefiltering is ever turned on for real work, this needs
a proper line-by-line translation of `FiltFiltM.m`'s ~218-line algorithm rather than relying on
scipy's built-in filtfilt variants.

## Design decisions

- **`fast_spectrum_batch` is one function for both auto- and cross-spectra**, matching the fact
  that `avgspec_sayan`/`avgspec_xy_sayan` in `mspec_fast.m` are mathematically identical (the
  auto version is just the cross version called with the same signal twice). Not reimplemented
  as two separate functions.
- **Complex-valued "floor" operation, auto-spectrum only**: MATLAB's `max(z, eps*maxz)` on a
  complex array compares by magnitude and keeps the actual complex element (not a real-part
  comparison). Replicated via `_complex_floor` in `fast_cross_spectrum.py`, but only for the
  auto-spectrum case (`x is y`) -- see "Known upstream bug" below for why the cross-spectrum case
  no longer applies it.
- **Trace ordering in the 2-D reshape** (`_reshape_to_traces_by_samples`): MATLAB's
  `reshape(A,[],N1)` on a column-major `(day,window,samples)` array produces traces ordered
  day-fastest-then-window. Reproduced in numpy via `transpose(2,1,0).reshape(n_samp,-1)` --
  chosen deliberately (not just "any" reshape that gives the right shape) so trace ordering
  matches exactly, which matters if per-trace diagnostics are ever compared 1:1 with MATLAB.
- **`S2`-before-`S1` argument order**: preserved exactly as in the original (only in the
  `'FastMspec'` branch) -- `sxy = fast_spectrum_batch(fmtse, s2, s1)`, i.e. X=S2, Y=S1. This is
  called out prominently in `crosscorr_mtc.py`'s docstring since it's easy to "fix" by mistake
  and it changes the conjugation/phase sense of the cross-spectrum.
- **One-sided output convention**: `fast_spectrum_batch` returns `floor(N/2)+1` bins (matching
  `avgspec_xy_sayan`'s convention), unlike `thomson_multitaper.FastMultitaper.spectral_estimate`
  which returns the full two-sided spectrum directly. This is intentional -- this pipeline does
  its own (now bug-fixed) one-sided-to-two-sided reflection afterward, matching the MATLAB
  source's structure, rather than reusing the single-signal class's different convention.

## Known upstream bug (fixed here, not reproduced)

The one-sided->two-sided reflection in `ccf_compute_crosscorr_mtc_Z.m`/`_T.m` checks the parity
of the already-one-sided array's length rather than the true window length `N`, giving a
wrong-length spectrum whenever `N mod 4` is 0 or 3. Fixed in `_reflect_onesided_to_twosided`
using `N`'s true parity. Full analysis: `../../verification/octave_verify_ccf_pipeline/README.md`. The one concrete
production config found in this codebase (`N=10801`, `N mod 4 == 1`) happens to avoid the bug,
so this wasn't caught by prior real-world use with that specific window length.

`avgspec_xy_sayan` (`functions/jSpectral/mspec_fast.m`, Sayan's cross-spectrum extension of
Karnik's `FastMultitaper.SpectralEstimate`) applies the same `max(z, eps*maxz)` positivity floor
that Karnik's original code uses for its (real, non-negative) auto-spectrum -- but `avgspec_xy_sayan`'s
`z` is `SXY`, a genuinely complex cross-spectrum with no non-negativity constraint. Lilly's
original (non-fast) `avgspec` never floors anything, for either the auto- or cross-spectrum case,
and Karnik's own code never produces a cross-spectrum at all, so there's no precedent in either
source for flooring a complex quantity this way. Because MATLAB's `max` on complex arguments
compares by magnitude, this silently overwrites any `SXY` bin near a coherence null -- exactly
where the phase carries the most physically meaningful information for Bessel-coherence/dispersion
work -- with a real value carrying the *global* maximum's phase instead of that bin's own. Fixed
in `fast_spectrum_batch` (`fast_cross_spectrum.py`) by only applying `_complex_floor` when `x is y`
(the auto-spectrum call convention used throughout `crosscorr_mtc.py`); the cross-spectrum case now
returns `(z0+z1)/K` unmodified. Not reproduced from the `.m` source. This is on the path exercised
by every real `'FastMspec'`-technique cross-correlation (`ccf_compute_crosscorr_mtc_T.m`/`_Z.m`),
so it affects `SXY` and therefore the coherence (`SXY/sqrt(SXX*SYY)`) used throughout Notebook 3.
**Expected test impact**: `tests/test_crosscorr_mtc.py::test_matches_octave_on_onesided_portion`
compares against a real Octave run of the unmodified (buggy) `.m` source -- it's skipped in
environments without the `octave_verify_ccf/synth_medium.mat` fixture, but wherever it does run
with that fixture, it will now legitimately fail on the FastMspec case's cross-spectrum comparison
near coherence nulls. That's expected given this is an intentional divergence, not a regression;
the fixture/test would need regenerating against a MATLAB source with the same fix applied (or the
test updated to tolerate this known, intentional difference) to pass again.

**Confirmed on real, large-scale data (2026-09-01, Notebook 5 revamp Stage 3)**: previously this
was only inferred from a small synthetic (N=32) Octave comparison. Recomputing FastMspec in Python
for the real SKRH-BAND pair (N=10801, 1605 traces) and comparing against Sayan's own precomputed
MATLAB `coh_sum` (`SKRH_BAND_fastmspec.mat`, generated by the unmodified, still-buggy `.m` source)
gives an overall one-sided relative L2 error of 3.06%, small in aggregate but concentrated exactly
where the floor-clamp signature predicts: binning by `|coh_sum|` (a proxy for distance from a
coherence null), relative error is 20.8% in the lowest-magnitude 10% of bins, 7.4% in the next 40%,
and 2.4% in the top 50% -- while the *absolute* error stays roughly constant (~1.05-1.12) across
all three bins. A roughly-constant absolute perturbation that only bites hard where the true signal
is small is precisely what a `max(z, eps*maxz)` floor produces, and rules out a genuine new
translation bug (which would more plausibly show a uniform or structurally different error
pattern). `coh_num` itself matches exactly (1605 both sides), confirming the trace-count/summation
logic is correct and the discrepancy is specifically in per-bin values, not bookkeeping.

## Real-data verification (2026-08-18, update)

Real SAC data became available (`data/test/raw_data/{SA53,SA58,MTAN,RUNG,...}`, `data/metadata/
{sta_list.txt,orientation.csv,...}`, plus a `bkup/` with a prior precomputed CCF result for the
SA53-SA58 pair). This closed out Phase 5's real-data-verification gap and Phase 6's end-to-end
gap. The SA53/SA58/MTAN/RUNG subset used here is available for download -- see the repo root
[`README.md`](../../README.md#getting-the-example-sac-data) ("Getting the example SAC data").

- `prepare_data.build_windows` on a real day-pair (SA58/SA53, `1998.289`, BHZ, winlength=4h)
  matches a real Octave run of the equivalent `.m` logic to ~1e-12 -- both the raw windowed
  samples and the resulting single-window coherency spectrum.
- **Full real-data pipeline** (`find_day_pairs` + `validate_pair` + `build_windows` +
  `dispatch.compute_crosscorr`'s plain-fft branch, across all 148 candidate day-pairs -> 70 valid
  days x 11 windows = 770 total, matching Octave's own count exactly) matches an independent full
  Octave run of the same real files to **3.9e-6 relative L2 error** -- the small residual is
  expected floating-point summation-order noise across 770 additions, not a translation defect.
- **A red herring worth recording**: the same comparison against the *archived* precomputed
  result in `bkup/results/.../SA58_SA53_f.mat` (from an experiment literally named
  `archive_expset4`) showed a large, structural mismatch (correlation ~0.02 across the bulk
  spectrum) -- initially alarming, but resolved by instead running Octave fresh on *today's*
  actual raw SAC files and finding it matches Python exactly. Conclusion: the archived reference
  was computed from a different (older/different) data snapshot than what's currently on disk,
  not evidence of a bug. Lesson: prefer a fresh Octave run on the exact same input files over an
  archived result whose provenance isn't certain, even when the archived result looks like
  "real ground truth" -- see `../../../verification/octave_verify_ccf_pipeline/real_data_test/` for the scripts used.
- The `IsMultiTaper=0/IsFTN=0/IsOBN=0/IsPrefilter=0/IsDetrend=0/IsTaper=0` config that produced
  this real reference exercises `dispatch.py`'s plain-fft branch, not the `IsMspec`/`FastMspec`
  path Phases 0-3 focused on -- so this is a valuable *independent* real-data check of a
  different code path than everything verified so far, not a re-check of the same one.

## What's genuinely unverified

- Only the Z-component function has been ported/tested (all three techniques now done).
- Verified against one synthetic fixture (`N=32`, random correlated data) -- not against real
  seismic data, and not against a range of `Wband`/`cutoff`/`epsilon` combinations.
- The `Iseval`/`Chunklen`/`zero_padding` performance-evaluation mode (from
  `ccf_synth_setup_perfrm_evl.m`) is not implemented -- only the non-eval path.
- The real-data check above (see "Real-data verification" section) only exercised the plain-fft
  dispatch branch (Phase 6) and `prepare_data.py` (Phase 5) end-to-end. The `IsMspec`/`FastMspec`
  path (Phases 0-3) is still only verified against synthetic data + Octave, not a real end-to-end
  SAC-to-output run -- worth doing now that real data and a known-good `prepare_data.py` are both
  available, e.g. by building `S1_data_mat`/`S2_data_mat` from real SA53/SA58 windows and feeding
  them into `compute_crosscorr_mtc_fastmspec` instead of synthetic random data.
- T-component (rotation via `orientation.csv`) hasn't been exercised against real N/E data yet,
  though the metadata file needed for it is now available.
