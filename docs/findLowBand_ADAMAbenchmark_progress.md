# findLowBand_ADAMAbenchmark: progress tracker

A new, separate notebook (sequenced after Stage 5 of the Notebook 5 revamp, not part of it) that
confirms or revises the bandwidth-selection theory (`NW_low`, `docs/stage5_bandwidth_theory.tex`)
against a **true ADAMA benchmark**: ADAMA's own real station pairs, with ADAMA's own real
`co`/`cf` dispersion curves as ground truth -- not Sayan's 380 (confirmed disjoint from ADAMA's
catalog, `docs/notebook5_revamp_progress.md`'s 2026-09-04/05 logs). Full scoping rationale for why
this is its own notebook, not folded into Stage 4.5 or Stage 5: see
`docs/notebook5_revamp_progress.md`'s 2026-09-05 "Scoping session" log entry.

## Goals (set 2026-09-05, user-directed)

1. Test different windowing strategies -- `N` (window length) should differ by station
   separation distance, per the already-derived Rayleigh/DFT-grid resolution bound
   (`NW_high(r) = n_bins/2`, `docs/round2_hypothesis_evaluation.tex` Section 5.1's "third
   constraint" subsection) and the window-length-vs-`NW_high` table already worked out there for
   Sayan's own quartiles. Not yet started for the ADAMA benchmark pairs specifically.
2. Verify ADAMA's full station-list coverage in whatever raw-data source we use (not just the XV
   subset) -- in progress, see Findings below.
3. Compare individual SAC files (terravibranium) against the compiled `ADAMA_gvib.h5`
   (repovibranium, `github.com/URseismology/ADAMA/tree/main/DataFiles`) -- done, see Findings.
4. Document the best data-access/processing strategy for this notebook -- this file.

## Findings (2026-09-05)

### Goal 2: station-list coverage on terravibranium

`ADAMA_stalist.csv` (62 networks, 1458 stations) vs.
`terravibranium:/RAID6/bluehiveBackup/Prj5_HarnomicRFTraces/2_Data/preprocessed_data/`:

- **Network-folder level: complete.** All 62 networks present as `Data{NET}/` folders, plus one
  extra (`IM`) not in the station list.
- **Station level, sampled not exhaustive**: `G` (6/6), `II` (8/8), `IU` (13/13), `YJ` (30/30)
  match exactly. `AF` does not: 42/47 present, missing `PWET`, `MONG`, `ZOMB`, `MZM`, `KIG` -- a
  real, quantified gap, not assumed away. Only 5 of 62 networks checked at station level so far;
  a full systematic check (all 62) is a concrete next step before treating this archive as
  complete, not yet done.
- The 6 mainland companion stations ADAMA's real Madagascar-touching pairs need (`MAPH`, `MOCU`,
  `MSGR`, `NAPU`, `SENA`, `TETE`) are present here, confirmed directly -- they are *not* present in
  Sayan's own bluehive archive (`PRJ_SPAC/data/test/raw_data/madagascar_data/DataXV/`, island
  stations only), which is why terravibranium's own ADAMA archive, not Sayan's, is the right
  source for this notebook.
- **Temporal overlap, spot-checked, not exhaustive**: real overlap confirmed between a handful of
  island/mainland station pairs (e.g. `MAPH` 2011.241-2013.241 vs. `BITY` 2012.217-2013.244;
  `DGOS` 2011.270-2013.236 fully contains `MAGY` 2012.074-2012.350). Only ~9 of the 34 XV stations
  checked; confirming overlap across all 141 island-touching pairs specifically is not yet done.

### Goal 3: individual SAC vs. `ADAMA_gvib.h5`

**Not directly comparable in practice, and the reason is concrete, not assumed.** `ADAMA_gvib.h5`
is 1.04 TB (`repovibranium:/volume1/web/ADAMA-D1/ADAMA_gvib.h5`). Tested whether a partial
download could still be opened for structural inspection (as worked earlier this session for the
much smaller `co`/`cf`/`ncfs_TT` files): pulled a 1.5 GB prefix, `h5py.File(...)` refused it
outright --
```
OSError: Unable to synchronously open file (truncated file: eof = 1885159424,
stored_eof = 1040173798832)
```
HDF5 needs the file's true trailing metadata to open at all; a prefix, however large relative to
the other files we've handled, is structurally useless here. A full download at this project's
observed relay throughput (~10 MB/s, repovibranium -> local -> bluehive) would take on the order
of a day, just to *inspect* the file -- not a reasonable cost for this comparison.

**Settled definitively by reading the actual writer code** (`DataReaderWriter/ppToHDF5.py`,
github.com/URseismology/ADAMA -- per direct user pointer), not left as inference. The entire
script:
```python
for path in glob.glob('/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/2_Data/preprocessed_data/Data*/*/*.sac'):
    stream = read(path)
    stream.write(namefile, 'H5', mode='a')
```
`ADAMA_gvib.h5` is built by globbing **exactly** `preprocessed_data/Data*/*/*.sac` -- the same
directory tree already found and verified on terravibranium -- reading each file with
`obspy.read()` and appending, no processing applied in this step. **`ADAMA_gvib.h5` and the
individual SAC files are the same data, byte-for-byte, just two different containers.** This also
means the earlier open question (raw counts vs. response-corrected) isn't resolved by this script
either way -- whatever the SAC files already are, `ppToHDF5.py` doesn't change it -- so that
question still stands, now clearly scoped to "whatever state the SAC files are already in," not
something the HDF5 packaging step could have altered.

**Useful related documentation, per direct user pointer**: `trichter/notebooks`'
`cross_correlation_okhotsk_coda.ipynb` (Tom Richter is `obspyh5`'s own author) demonstrates the
proper `obspyh5` API for this kind of file: `obspyh5.set_index(...)` configures the internal
group-naming scheme before writing (confirms the `index` file-attribute already found in
`ADAMAraw_co_love.h5` is a standard, documented mechanism, not an ADAMA-specific hack), and
`from obspyh5 import iterh5; iterh5(path)` **streams traces one at a time** rather than requiring
the whole file in memory. Doesn't change the recommendation below (the real bottleneck is the 1 TB
*download*, not in-memory loading once downloaded) but is the right tool if `ADAMA_gvib.h5` is
ever used directly in a later stage.

**One open question still flagged, not resolved**: whether terravibranium's SAC files are raw
instrument counts or already response-corrected. `ccf_prepare_data_T_mdg.m` references a
`PZpath` (pole-zero response files) as a separate input, suggesting Sayan's own pipeline expects
to do its *own* response removal -- i.e. the SAC files are likely still raw counts despite living
under a folder named `preprocessed_data`. Not directly confirmed (would need to inspect an actual
SAC header, e.g. `IDEP`, or check for a companion RESP/PZ archive for the mainland stations) --
flagged for whoever starts building the actual SAC-to-mat pipeline, not assumed either way.

## Recommendation (Goal 4) -- revised, `ADAMA_gvib.h5` is the decision

My first pass here (below, struck through in substance not in this text) recommended the
individual SAC files over `ADAMA_gvib.h5`, weighted mainly on download size/practicality.
**Overridden by direct user decision**, on grounds I'd underweighted: `ADAMA_gvib.h5` is a single
file (simpler to reference than a directory tree of many thousands of small SAC files), supports
genuine parallel access (multiple readers can open one HDF5 file concurrently -- important for a
SLURM array of NW-sweep workers, and avoids the small-file metadata overhead a parallel filesystem
like bluehive's scratch pays for opening/listing millions of individual files), and needs no
low-level directory-hierarchy parsing (`glob`-ing a big nested tree) at read time. These are real,
durable advantages for a notebook that will read this data repeatedly, not just once -- the
one-time cost of acquiring the file is worth paying for that.

**Status, corrected**: `ADAMA_gvib.h5` was already local on bluehive the whole time -- a first
`find` was killed mid-run before completing (per the "stop crawling" correction) and read as
"nothing found," which was wrong; the same `find`, left to finish, located it at
`/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5`. Confirmed via exact
byte size (`1,040,173,798,832` -- identical to repovibranium's copy): this is the complete file,
not a partial or older version. A ~1.04 TB cross-network transfer was started before this was
known and has been killed, and the partial download deleted -- **no transfer is needed at all**;
this notebook can read the file directly from where it already sits. Lesson worth stating
plainly: killing an in-progress investigation early can produce a confident-sounding but wrong
negative result -- the earlier "confirmed absent" claim should have waited for the check to
actually finish, or been stated as "not found in the time available" instead.

**Bonus finds in the same directory, inspected (cheap now -- local file, full random access, no
transfer needed)**: `ADAMA_gvib_failed.h5` (593 GB) and `ADAMA_gvib_test.h5` (550 KB, likely a
development-time sample, not further inspected). `_failed` is a real `obspyh5` waveform archive
too, same indexing scheme, **1372 stations in the main file (matches ADAMA1's paper's "~1372
stations" exactly), 857 in `_failed` -- and every one of those 857 is also present in the main
file** (not a disjoint missing-stations list). So `_failed` isn't "stations we don't have" -- it
holds specific day/window traces that failed some later processing step for stations that
otherwise have good data in the main archive. All 33 of our XV stations appear in both files
identically, meaning: none of our relevant stations are entirely missing or entirely failed, only
possibly some individual days within them. This doesn't fully replace the function-based check
above (a station having *some* good data doesn't guarantee the *specific* days needed for a given
pair overlap correctly), but it rules out the worst-case failure mode (a station silently absent
from the archive) directly and cheaply.

**A real risk in `ADAMA_gvib.h5`, worth stating plainly, not assumed away**: `ppToHDF5.py`'s own
error handling is `try: ... except: print(path)` -- it silently swallows failures, so some SAC
files could in principle be missing from `gvib.h5` without any retained record. **The right check
for this is not an inventory audit** (comparing every SAC path against every `gvib.h5` entry would
itself mean crawling the SAC tree -- exactly what's already been ruled out as unproductive).
**Correction, per direct user guidance**: verify by *function*, not by *inventory* -- pull a known
pair's data back out of `gvib.h5` (e.g. `AF.SKRH`-`XV.BAND`, already validated in Stage 3 against
Sayan's own precomputed MATLAB `coh_sum` -- `verification/skrh_band_real_data/`), run it through
the same windowing/cross-correlation pipeline, and confirm the result matches that already-known-
good baseline. If it recovers correctly for known pairs, `gvib.h5` is trustworthy for this
purpose regardless of whether every single SAC file in the tree happens to have made it in --
next concrete task after the transfer lands.

**After that verification, the actual deliverable (per direct user request, explicit that prior
documentation here wasn't sufficient)**: a properly designed, well-documented Python library --
not ad-hoc scratch scripts -- for loading *pair-matched* data out of `gvib.h5` (which is indexed
per single-station-channel-day, not per-pair, so the library must find the actual overlapping-day
set for two given stations itself), with correct windowing and overlap (parameterized, not
hardcoded to Sayan's fixed 3-hour convention -- Goal 1 needs `N` to vary by distance, so the
library's windowing must be a real parameter, not a constant), and designed for safe concurrent
access (multiple read-only `h5py` handles on the same file from parallel workers) so a SLURM-array
NW-sweep can actually use it. Matches the documentation depth already set by
`python/dispcurve_pick/hybrid_reference_curve.py` -- module docstring covering method and
provenance, a companion README covering data/validation -- as the bar to clear, explicitly, not
to fall short of again.

**Before running the actual SAC-to-mat/pair-loading logic**, beyond the 1:1 verification above:
1. Confirm the raw-counts-vs-response-corrected question (cheap: inspect one SAC header).
2. Extend the temporal-overlap check from a spot sample to all 141 island-touching pairs.
3. Decide the windowing/overlap convention this library exposes as parameters -- reusing
   `ccf_prepare_data_T_mdg.m`'s exact convention as the *default*, matching this project's
   established practice of verifying a Python reimplementation against the original MATLAB/Octave
   output before trusting it (`python/ccf_pipeline/NOTES.md`'s own pattern), while allowing Goal 1's
   distance-dependent `N` to override it.

## Goal 1: windowing strategy (design note, not yet implemented)

`NW_high(r) = n_bins/2` (the Rayleigh/DFT-grid bound) ties directly to window length `N`: more
bins means either a longer window (larger `N` at fixed sample rate) or a higher sample rate at
fixed duration. `docs/round2_hypothesis_evaluation.tex` already has a concrete worked table
(window-length-vs-`NW_high` by quartile, Section 5.1's "concrete, quantitative fix" subsection) as
a starting point for how `N` should scale with distance `r` -- this notebook's job is to test that
scaling directly on ADAMA's own real pairs (which span a real, if narrower, distance range:
724-2220 km for the 141 island-touching pairs) rather than assume the existing table transfers
unchanged. Not yet started.
