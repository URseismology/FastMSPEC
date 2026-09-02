# `dispcurve_pick_batch`: the 380-pair x 4-technique bluehive batch pipeline

[← Back to repo README](../../README.md) | See also: [python/dispcurve_pick/NOTES.md](../dispcurve_pick/NOTES.md) | [verification/skrh_band_real_data/README.md](../../verification/skrh_band_real_data/README.md) (the small-scale validation this generalizes) | [docs/notebook5_revamp_progress.md](../../docs/notebook5_revamp_progress.md) (Stage 4 log)

Recomputes all 380 Madagascar station pairs across all 4 `ccf_pipeline` techniques
(single-taper, FastMspec, Mspec, MspecBestK), picks a phase-velocity curve for each via the
instrumented `dispcurve_pick` picker scanned across the `SDISPL.ASC` +/- 0.8 km/s template
corridor, and cross-validates FastMspec/single-taper against Sayan Swar's own precomputed MATLAB
cross-spectra where available (363/380 pairs). Everything here was validated at small scale first
(`verification/skrh_band_real_data/`) before being generalized to a full batch.

## Layout

- `manifest.py` -- reads the 380-pair catalog, builds the full 1520-work-unit list
  (`(pair, technique)`), **technique-outer ordering** so each technique occupies one contiguous
  index range (`technique_index_ranges()`) -- Mspec costs >10x every other technique (Stage 3
  timing), so a single uniform SLURM resource request would be wrong for 3/4 of the work units;
  submitting one array job per technique with technique-sized resources needs contiguous ranges.
- `work_unit.py` -- the core `process(pair, technique, ref_curve_path) -> WorkUnitResult`: loads
  the matched-data `.mat`, computes the cross-spectrum with the right preprocessing for that
  technique (single-taper needs detrend+5%-cosine-taper *first*; the other three use the
  un-preprocessed matched data directly -- see `verification/skrh_band_real_data/README.md` for
  why), optionally cross-validates against a precomputed MATLAB reference, scans the template
  corridor through the instrumented picker, and keeps the best-scoring template. Never raises --
  any exception is caught and recorded in the result's `error` field, so one bad pair/technique
  can't take down a batch of 1520.
- `run_plain.py` / `run_multiprocessing.py` -- two drivers over the same `work_unit.process()`
  core (see "Plain vs. multiprocessing" below). Each work unit's result is written as its own
  small JSON file (`<pair_id>__<technique>.json`) in the results directory, not appended to a
  shared CSV -- many concurrent SLURM tasks writing to one file is a race condition waiting to
  happen. Both are idempotent/resumable: an existing result file is skipped, not recomputed.
- `aggregate.py` -- combines whatever result JSON files exist into one `manifest.csv`. Safe to run
  on a partial batch (a partial, valid manifest comes out); re-run any time to pick up more.
- `submit_plain.sbatch` / `submit_mp.sbatch` -- SLURM array job templates, one submission per
  technique (see "Per-technique array ranges" below).

## Plain vs. multiprocessing

Per direct guidance during Stage 3 (observed directly: running all 4 techniques serially in one
process is slow and memory-heavy, no benefit from serializing independent computations), the plan
is to build both and decide from measured cost, not assumption:
- `run_plain.py`: one work unit per invocation -- maps 1:1 onto a SLURM array task
  (`SLURM_ARRAY_TASK_ID` selects the work unit). Simple, maximally fault-isolated (one task's
  failure/preemption never affects another), but pays fresh Python/numpy/scipy/obspy import
  overhead per work unit and can leave a node's other cores idle if only one array task lands on
  it at a time.
- `run_multiprocessing.py`: one invocation processes a *slice* of work units, fanned out across a
  `multiprocessing.Pool` sized to `--cpus-per-task`. Maps onto one SLURM array task claiming a
  full node, using all its cores via one long-lived interpreter per worker (import cost paid once
  per worker, not once per work unit).

Both exist specifically so the choice between them (and the right chunk size / worker count) can
be made from real measured throughput, per direct guidance -- see
`docs/notebook5_revamp_progress.md`'s Stage 4 log for the actual local-machine comparison run
before this was ever deployed to bluehive, and for whichever driver the real bluehive run ended
up using.

## Per-technique array ranges

From `manifest.technique_index_ranges()` (380 pairs each): single-taper `[0, 380)`, FastMspec
`[380, 760)`, Mspec `[760, 1140)`, MspecBestK `[1140, 1520)`. Submit one array job per technique
with `--array` restricted to that range and `--time`/`--mem` sized for that technique specifically
(see the comment blocks at the top of each `.sbatch` file for concrete example invocations) --
Mspec alone needs roughly the time/memory budget of all three other techniques combined, per
Stage 3's timing pilot (~45 min cross-spectrum vs. ~2.5 min for FastMspec/MspecBestK, ~2.5s for
single-taper).

## Deploying to bluehive

**Deployed location differs from the original plan**: `PRJ_SPAC/codes/prod/` (inside Sayan's own
`Sayan_Swar_WS`) turned out to be read-only to this account (`drwxr-x---`, group has no write) --
discovered live, not assumed. New work instead lives at a fresh top-level directory,
`/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/` (that account owns `/scratch/tolugboj_lab/`
itself, confirmed writable) -- still clearly labeled, still nothing else touched, just not nested
inside another lab member's personal workspace.

**Environment differs from the original plan too**: a plain `venv` with `pip install` hit a wall
of genuine toolchain fragility on the login node -- the default `python3/3.11.0` module's SSL
module doesn't work at all (no OpenSSL support compiled in, every HTTPS pip request silently
"finds no versions"); `python3/3.11.10` fixes that, but its `python -m venv` needs
`LD_LIBRARY_PATH` from that same module load, which (like `$BASE` env vars generally) doesn't
persist across separate `ssh` invocations -- every command touching the venv needs `module load
python3/3.11.10` in the *same* command; building numpy/scipy from source needs a newer gcc than
the system default `cc` (`module load` provides one, but meson invokes `cc` specifically --
`CC=gcc CXX=g++` fixes that); and `obspy`'s `pyproj` dependency needs PROJ >= 9.4.0, newer than
any `proj` module available (max 8.1.1, itself broken by an unrelated OpenSSL library-version
mismatch). None of this is a dead end, but stacking workaround on workaround for pure-Python-ish
packages was the wrong tool for the job here. Switched to a **fresh, dedicated conda environment**
(`fastmspec_batch`, in the same shared Anaconda install other lab members use --
`/scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/` -- but a new env, not Sayan's own
broken `Seislib` one, which stays untouched) instead: `conda create -n fastmspec_batch -c
conda-forge python=3.11 numpy scipy pandas obspy` installed everything as prebuilt binaries, no
compilation, no toolchain issues at all. Still "a fresh, self-contained environment, not the
shared broken one" in spirit -- just conda as the packaging mechanism instead of venv, since venv
turned out to demand exactly the compilation this login node can't reliably do.

**`seislib` itself is not installed on bluehive at all, and doesn't need to be** -- see
`python/dispcurve_pick/NOTES.md` "Why the `seislib` package dependency was removed": the picking
path was decoupled from the full `seislib` package specifically because `pip install seislib`
also fails to build here (a genuinely broken, unrelated Cython extension in `seislib.tomography`,
independent of the numpy/pyproj issues above) -- vendoring the ~6 small functions/classes actually
needed sidesteps this entirely, on top of being a good idea regardless of platform.

Steps, in order:
1. ~~Vendor packages into a venv~~ -- see above; use the `fastmspec_batch` conda env instead:
   `conda create -y -n fastmspec_batch -c conda-forge python=3.11 numpy scipy pandas obspy`.
2. Copy this repo's `python/ccf_pipeline/`, `python/thomson_multitaper/`, `python/dispcurve_pick/`,
   `python/dispcurve_pick_batch/`, plus `data/reference/SDISPL.ASC` and the pulled
   `madagascar_stn_conn_ccflist.csv`, into `/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/`.
3. Validate on `debug` (12 nodes, 1-hour cap) with a tiny 1-2-work-unit array first -- confirms
   the SLURM mechanics and paths work before committing real partition time.
4. Run the real batch on `preempt` (see the plan's `sinfo` investigation for why), with `urseismo`
   as fallback/supplement once its current jobs free up.

The `.sbatch` templates' `source env/venv/bin/activate` line is stale (written before this
environment change) -- update to `source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/
2021.05/etc/profile.d/conda.sh && conda activate fastmspec_batch` before submitting anything for
real; not yet fixed in the committed templates as of this note (see
`docs/notebook5_revamp_progress.md`'s Stage 4 log for current status).

## Cross-validation coverage

Only `FastMspec` and `single-taper` have precomputed MATLAB references (363/380 pairs each), at
`PRJ_SPAC/results/test/love/madagascar/{fastmspec,firstorder}/ccf/window3hr/fullStack/ccfTT/
<STA1>/<STA1>_<STA2>_f.mat` (confirmed by direct listing -- 363 files each). `Mspec`/`MspecBestK`
have no MATLAB reference anywhere in this project; their results are new, Python-only findings,
reported as such in the eventual manifest/notebook, not silently treated as equally verified.
