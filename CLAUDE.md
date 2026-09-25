# CLAUDE.md — FastMSPEC (agent entry point)

Read this file first, then [`docs/STATE.md`](docs/STATE.md) (what is done / open, per work stream), then only the
one stream document you need. [`docs/INDEX.md`](docs/INDEX.md) lists every document and whether it is current,
historical or superseded. Do **not** start by reading `docs/notebook5_revamp_progress.md` end to end: it is a
2,285-line dated log; use its "Status summary" section (top) and search the log by date when you need a reason.

## What this project is
Python translation and validation of FastMSPEC (fast multitaper cross-spectra, Karnik–Romberg–Davenport 2022) for
ambient-noise cross-correlation, built for Sayan Swar's Madagascar / ADAMA work. Current program goal (Tolu,
2026-09-25): validate FastMSPEC **at scale** — ~2,000 stations worldwide, global Love-wave noise correlations, global
dispersion maps compared with Ekström (GDM52), demonstrating FastMSPEC's improvement. See
[`docs/global_validation_readiness.md`](docs/global_validation_readiness.md).

## Layout (one line each)
| path | contents |
|---|---|
| `python/thomson_multitaper/` | multitaper library (`FastMultitaper`, `Multitaper`, adaptive); verified against Octave |
| `python/ccf_pipeline/` | translation of the CCF multitaper cross-correlation pipeline (techniques: single-taper, Mspec, MspecBestK, FastMspec) |
| `python/dispcurve_pick/` | vendored + instrumented `seislib` picker, hybrid ADAMA+GDM52 reference curve, `gvib_loader.py` (reads `ADAMA_gvib.h5`) |
| `python/dispcurve_pick_batch/` | bluehive batch pipeline (380 pairs × 4 techniques; Round 2 bandwidth sweep) |
| `notebooks/` | six teaching/validation notebooks; sources are `notebooks/_lib/build_nb*.py`, `.ipynb` files are generated and committed with outputs |
| `verification/` | every check that a translation or claim is right (Octave comparisons, real-data checks, Notebook 3 precompute drivers) |
| `docs/` | trackers, theory (`.tex`/`.pdf`), reports, references — see `docs/INDEX.md` |
| `data/` | mostly gitignored; small committed results under `data/results/`, reference curves under `data/reference/` |
| `legacy/` | original MATLAB source and the repository's pre-translation contents (read-only history) |

## Rules that are easy to get wrong
1. **Notebooks are generated.** Edit `notebooks/_lib/build_nb<N>.py` (and `nb3_helpers.py` / `nb4_helpers.py`), regenerate, and
   **execute for real** with `nbclient` (`NotebookClient(nb, timeout=1800, kernel_name="python3")`); never hand-edit
   `.ipynb` JSON except as a documented surgical fix, and never commit a zero-output notebook as "current".
2. **Real data, not toys.** Findings must be shown on real station pairs (best/hardest example, "show before telling").
   Do not claim more than the numbers support; state failures with exact numbers.
3. **Use full words** in prose (write "peak memory", "received level", "transmission loss"; avoid unexplained abbreviations).
4. **Round 1 / Round 2 results used the pre-Stage-4.5 picker** (single generic `SDISPL.ASC`, uniform corridor,
   `horizontal_polarization=False`). They are the historical record, not current-pipeline output.
5. **Love vs Rayleigh:** Love waves need `horizontal_polarization=True` (J0−J2 model); vertical-component (Rayleigh) uses J0. Only Love has
   been validated against ADAMA in the picker; the Rayleigh path is unvalidated.
6. **Resolution ceiling:** `NW_high = N·c_min / (4R)` (window length N in samples, minimum velocity `c_min`, station separation R),
   with a numerical floor near NW ≈ 3 (K = 0). Floor `c_min` at 3.0 km/s when the hybrid curve gives an unphysical minimum.
7. **Data sources:** `ADAMA_gvib.h5` (1 TB, on bluehive only) is the decided data source for ADAMA benchmarks; SAC files are
   already instrument-response-corrected — do not reopen either question. Sayan's matched-data `.mat` files (380 pairs) are a
   different, MATLAB-path dataset.
8. **Do not edit generated or vendored files casually:** `python/dispcurve_pick/_vendored_seislib_*.py` stays byte-identical to upstream
   except the documented instrumentation (see `python/dispcurve_pick/NOTES.md`).

## Environments
* **Local:** system `python3` (numpy, scipy, obspy, pytest; notebooks also need `nbformat nbclient ipykernel`, a registered `python3` kernel, and
  `seislib` for `python/dispcurve_pick/tests/test_matches_upstream.py`). Quick check:
  `cd python && python3 -m pytest -q --ignore=dispcurve_pick/tests/test_matches_upstream.py` → 20 passed, 12 skipped (2026-09-25).
* **Bluehive** (SLURM; account `tolugboj_lab`, partition `urseismo`, or `standard` with `--qos=standard`): batch root
  `/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch`, conda env `fastmspec_batch`. The controller times out intermittently: before resubmitting after a
  failed `sbatch`, check `sacct --name=<job>`; poll `sacct` states rather than `squeue`.
* **Local data** (gitignored; re-pull per [`data/README.md`](data/README.md) and `data/reference/hybrid_curve_README.md`): SKRH–BAND matched-data
  `.mat`, ADAMA maps + GDM52, the SAC example data.

## Git conventions
Branch `notebook5-phase-velocity-revamp` is the working branch (`main` stops at `c73eecd`, 2026-09-04, Round 2 complete — 34 commits behind).
Commit early, meaningful messages; commit messages end with `Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>`.
**Ask before pushing** this repository. Tag `notebook5-v1-event-scanning` preserves the pre-renumbering notebooks.

## Adjacent projects (outside this repository, local paths on Tolu's machine)
* `~/claude-sandbox/projects/wavenet_xd_pair_test/` — smoke test of the WaveNet-EpicAI acquisition pipeline
  (`URseismology/wavenet-epicAI`) on XD.MTAN/XD.RUNG plus verified fixes. It gates the global run; living handoff is its `STATE.md`.
  Do not touch `wavenet_ncf_production` or `wavenet_ncf_canary` on bluehive.

## When you finish a piece of work
Update the relevant row in `docs/STATE.md` (status + date + commit), add a dated entry to the stream's own tracker if it has one,
and run `python3 docs/check_links.py` so the documentation stays navigable.
