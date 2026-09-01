# SKRH-BAND real-data verification

[← Back to repo README](../../README.md) | See also: [python/ccf_pipeline/NOTES.md](../../python/ccf_pipeline/NOTES.md) | [python/dispcurve_pick/NOTES.md](../../python/dispcurve_pick/NOTES.md) | [docs/notebook5_revamp_progress.md](../../docs/notebook5_revamp_progress.md) (Stage 3 log)

The first real, large-scale (N=10801, 1605 traces) end-to-end check of both `ccf_pipeline`'s
FastMspec/single-taper cross-spectrum computation and `dispcurve_pick`'s instrumented picker --
previously each was only verified against small synthetic fixtures (`ccf_pipeline`) or the
installed `seislib` package on synthetic input (`dispcurve_pick`). SKRH-BAND (290.3 km) is also
the exact pair Sayan Swar's own executed notebook (`phasevel_compute_slide13and14.ipynb`, on the
lab NAS) picked a dispersion curve from, giving this a real known-outcome to reproduce rather than
only an internal-consistency check.

## Getting the data (not committed -- matches this repo's established practice for bulk/derived data)

Three files, none committed here:
```bash
# Raw matched windowed data (267MB) -- from bluehive
ssh bluehive "cat /scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/test/processed_data/love/madagascar/AFSKRH_XVBAND_win_3_all_matched_data.mat" > AFSKRH_XVBAND_win_3_all_matched_data.mat

# Precomputed MATLAB cross-spectra (164KB each) -- from the lab NAS
ssh repovibranium "cat /volume1/web/FastMSPEC_data/madagascar_data/pre_computed_files/SKRH_BAND_fastmspec.mat" > SKRH_BAND_fastmspec.mat
ssh repovibranium "cat /volume1/web/FastMSPEC_data/madagascar_data/pre_computed_files/SKRH_BAND_firstorder.mat" > SKRH_BAND_firstorder.mat
```
Both hosts' SFTP/scp subsystems are restricted to certain shared folders; plain `ssh ... cat`
works around it (same pattern already used for `data/reference/SDISPL.ASC` and
`docs/references/`'s papers). Then: `python3 validate_skrh_band.py`.

## What it checks, and what it found

**Cross-spectrum cross-validation** (Python recompute vs. Sayan's own precomputed MATLAB
`coh_sum`, same pair):
- FastMspec: ~3.06% relative L2 error, concentrated almost entirely near coherence nulls (20.8%
  in the lowest-magnitude 10% of bins vs. 2.4% in the top 50%, while absolute error stays roughly
  constant across bins) -- the exact signature of the already-documented, intentionally-not-
  reproduced MATLAB complex-floor bug (`ccf_pipeline/NOTES.md`'s "Known upstream bug" section).
  Confirms that bug's existence and behavior on real data for the first time; not a regression.
- single-taper: must apply detrend + 5% cosine taper *before* the plain-FFT coherency (this
  project's established "5% Cosine Single-Taper" technique, e.g. Notebook 3 Section 2) -- not raw,
  unprocessed data straight into FFT. Once applied: ~4.8e-9 relative error, machine precision.
  Omitting that preprocessing gives a misleadingly large ~46.5% mismatch that looks like a bug in
  the pipeline but is actually a missing preprocessing step in the *caller*.

**Dispersion-curve picking** (reproducing `phasevel_compute_slide13and14.ipynb`'s own real,
executed result -- checked directly from that notebook's saved outputs, not assumed): FastMspec
converges to a clean picked curve (~0.02-0.115 Hz, ~3.4-3.85 km/s); single-taper does not, even
after extra `savgol_filter` smoothing at a narrower band. Reproduced exactly with the vendored+
instrumented picker, plus a fair symmetric comparison (both raw, same band) Sayan's own notebook
never tried -- single-taper fails even harder unsmoothed (2066 candidate crossings, only 15
accepted picks, 0% frequency coverage). One design-relevant finding: `bad_quality_fraction` alone
is a weak discriminator across these three runs (~0.67-0.72 throughout) -- `n_accepted_picks` and
`freq_coverage_fraction` are what actually separate a converging pick from a failing one.

## Timing (provisional, this local machine -- see `docs/notebook5_revamp_progress.md` Stage 3 for
the full numbers and Stage 4's partition/parallelism design built on them)

Mspec dominates by >10x over every other technique (~45 min vs. ~2.5 min for FastMspec/MspecBestK,
~2.5s for single-taper) -- needs its own, more generous task budget in any batch design, not
uniform sizing across techniques.
