# FastMSPEC

Python translation of a MATLAB ambient-noise cross-correlation (CCF) pipeline built on fast
multitaper spectral estimation, based on Santhosh Karnik's "Fast Slepian Transform" method
(S. Karnik, Z. Zhu, M. B. Wakin, J. K. Romberg, and M. A. Davenport, "The fast Slepian
transform," *Appl. Comput. Harmon. Anal.* 46(3):624-652, 2019,
[doi:10.1016/j.acha.2017.07.005](https://doi.org/10.1016/j.acha.2017.07.005)). See
[`notebooks/README.md`](notebooks/README.md#references) for the full reference list, including
the follow-up multitaper paper and the ambient-noise application this repo was built for.

This repository was repurposed and renamed from `Para_CCF` (2021), the earliest parallel-CCF
MATLAB codebase in this project's lineage — see [`legacy/`](legacy/) for how the two connect.

## What's here

| Directory | Contents |
|---|---|
| [`python/`](python/) | The Python translation: [`thomson_multitaper/`](python/thomson_multitaper/) (DPSS/multitaper spectral estimation library, notes in [`python/NOTES.md`](python/NOTES.md)) and [`ccf_pipeline/`](python/ccf_pipeline/) (the ambient-noise CCF pipeline that calls it, notes in [`python/ccf_pipeline/NOTES.md`](python/ccf_pipeline/NOTES.md)). |
| [`verification/`](verification/) | The GNU Octave environments used to verify the Python translation against the real, unmodified MATLAB source — not synthetic reimplementation, the *actual* `.m` files, run and diffed against Python output. See [`octave_verify_multitaper/README.md`](verification/octave_verify_multitaper/README.md) and [`octave_verify_ccf_pipeline/README.md`](verification/octave_verify_ccf_pipeline/README.md). |
| [`legacy/matlab_source/`](legacy/matlab_source/) | The original MATLAB source this project translates from and verifies against, for side-by-side reading. See its [`README.md`](legacy/matlab_source/README.md). |
| [`legacy/para_ccf_original/`](legacy/para_ccf_original/) | This repo's original contents before the rename — the 2021-era predecessor scripts. See [`NOTE.md`](legacy/para_ccf_original/NOTE.md) and the original [`README.md`](legacy/para_ccf_original/README.md). |
| [`python/dispcurve_pick/`](python/dispcurve_pick/) | A vendored, instrumented copy of `seislib`'s dispersion-curve picker, plus the hybrid ADAMA+GDM52 per-pair reference-curve library and the gvib.h5 pair-matched data loader. |
| [`python/dispcurve_pick_batch/`](python/dispcurve_pick_batch/) | The bluehive batch pipeline that runs the picker across the full 380-pair dataset (multiple techniques, an NW-bandwidth sweep). |
| [`docs/`](docs/) | The technical plan for the CCF pipeline translation ([`plan_ccf_mtc_translation.md`](docs/plan_ccf_mtc_translation.md)), and — the current active work stream — the Notebook 5 revamp ([`notebook5_revamp_progress.md`](docs/notebook5_revamp_progress.md), **read its "Status summary" section first**) and its planned follow-on ([`findLowBand_ADAMAbenchmark_progress.md`](docs/findLowBand_ADAMAbenchmark_progress.md)). |
| [`notebooks/`](notebooks/) | Theory-to-application documentation, six planned Jupyter notebooks: why multitaper spectral estimation and why "Fast" (reproducing Karnik et al.'s own paper figures), why this matters for ambient-noise cross-correlation, the pipeline applied to real data, phase-velocity dispersion-curve picking at scale (mid-rebuild), bandwidth-selection theory vs. a real ADAMA benchmark (not yet started), and a scoped future-work roadmap for coda-correlation. See [`notebooks/README.md`](notebooks/README.md). |

## Status, in one paragraph

Both packages are translated and verified against the real MATLAB source using GNU Octave —
not just independent reference implementations, but the actual `.m` files, run directly and
diffed against Python output line-by-line. Three real bugs were found in the upstream MATLAB
along the way (two in `firstNlambdaDPSS.m`, one in a spectrum-reflection step shared by the CCF
functions) and fixed, with before/after Octave runs confirming each fix. The full pipeline —
real SAC file loading through to a computed cross-correlation — has been run end-to-end on real
seismic data (station pair SA53/SA58) and matches an independent Octave run of the same real
files to a relative error of 3.9e-6. Open items and exact scope boundaries are listed in each
package's `NOTES.md` — nothing here is claimed more thoroughly verified than it actually is.

## Current work: Notebook 5 revamp (phase-velocity picking)

This CCF pipeline now feeds a second, actively-developed work stream: judging pick quality by
tracking a full phase-velocity dispersion curve (an instrumented, vendored copy of `seislib`'s
picker), not just scanning raw zero-crossings. **Source of truth for current status is
[`docs/notebook5_revamp_progress.md`](docs/notebook5_revamp_progress.md)'s own "Status summary"
section** (read that, not this paragraph, for anything beyond the one-line summary here) — as of
2026-09-08: the picker is vendored and validated, the full 380-pair dataset has been run and
analyzed on bluehive (Round 1 + a 300-point bandwidth sweep, `docs/round2_hypothesis_evaluation.pdf`),
and a per-pair hybrid reference curve + corridor-search fix (Stage 4.5) is built, validated, and
merged into the production pipeline. Next: rebuilding the dispersion-curve notebook itself (Stage 5,
`notebooks/04_dispersion_curve_picking.ipynb` — notebooks were renumbered 2026-09-08, see
`notebooks/README.md`), then a follow-on notebook (`05_bandwidth_selection.ipynb`) benchmarking
against ADAMA's own real station pairs
([`docs/findLowBand_ADAMAbenchmark_progress.md`](docs/findLowBand_ADAMAbenchmark_progress.md)).

## Where to start reading

1. [`python/NOTES.md`](python/NOTES.md) and [`verification/octave_verify_multitaper/README.md`](verification/octave_verify_multitaper/README.md)
   — the multitaper library and how it was verified.
2. [`docs/plan_ccf_mtc_translation.md`](docs/plan_ccf_mtc_translation.md) — why the CCF pipeline translation targets
   `ccf_compute_crosscorr_mtc_Z.m`/`_T.m` specifically (not the originally-assumed
   `a1_ccf_ambnoise_*.m` scripts), and the phased plan that followed from that finding.
3. [`python/ccf_pipeline/NOTES.md`](python/ccf_pipeline/NOTES.md) and [`verification/octave_verify_ccf_pipeline/README.md`](verification/octave_verify_ccf_pipeline/README.md)
   — the CCF pipeline translation itself, phase by phase, including the real-data verification.
4. [`legacy/matlab_source/README.md`](legacy/matlab_source/README.md) — a map of which original `.m` file backs which Python
   module, for line-by-line comparison.
5. [`legacy/para_ccf_original/NOTE.md`](legacy/para_ccf_original/NOTE.md) — how this repo's original 2021 contents connect to the rest.
6. [`notebooks/`](notebooks/) — if you want the *why*, not just the *what*: theory, motivation, and real-data
   application, worked through with real figures rather than summarized.

## Getting the example SAC data

The real seismic data used in the tests, examples, and real-data verification (station pairs
SA53/SA58 and MTAN/RUNG, ~1.2 GB) is hosted on the lab's Synology NAS rather than committed to
this repo or a public GitHub release, since it's real network waveform data rather than something
appropriate to redistribute broadly:

```bash
curl -L https://repovibranium.synology.me/FastMSPEC_data/raw_data.tar.gz | tar xz -C data/
```

This unpacks into `data/raw_data/` and `data/metadata/`, matching the paths the tests and
examples expect. Requires access to the lab network/VPN. SHA256 of the tarball:
`4d5e4b3124d5d1325f57618d23d38755ba2c41ffb0a0c1711d3d8d0dc9d3ef46`.

## Running the tests

```bash
cd python
python3 -m pip install -r requirements.txt
PYTHONPATH=. python3 -m pytest tests/ ccf_pipeline/tests/ -v
```

Some `ccf_pipeline` tests compare against `.mat` fixtures generated by the Octave scripts in
[`verification/octave_verify_ccf_pipeline/`](verification/octave_verify_ccf_pipeline/) and will skip (not fail) if those fixtures aren't
present — see that directory's [`README.md`](verification/octave_verify_ccf_pipeline/README.md) to regenerate them (requires GNU Octave with the
`signal` package installed).
