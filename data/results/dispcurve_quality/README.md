# Stage 4 Round 1 results — dispersion-curve pick quality across the full dataset

`manifest.csv` (1520 rows) + `results/*.json` (one file per row) — the complete, aggregated
output of the Notebook 5 revamp's Stage 4 batch: all 380 station pairs x 4 techniques
(single-taper, FastMspec, Mspec, MspecBestK), run on bluehive.

**Provenance**: computed by `python/dispcurve_pick_batch/` (`work_unit.py`'s `process()`),
submitted via `submit_plain.sbatch`/`submit_mp.sbatch` across jobs `31344560`-`31347891` (see
`docs/notebook5_revamp_progress.md`'s Stage 4 log for the full, dated account of every bug found
and fixed along the way — 2D-shape collapse, FastMspec runtime variance, Mspec's systemic
mis-sizing, MspecBestK OOM on the largest files, and 7 leftover pre-fix memory errors caught only
by the final aggregation pass). Aggregated via `python3 -m dispcurve_pick_batch.aggregate`, pulled
from `/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/` via `scp` (888KB total — small enough to
commit directly, unlike the raw ~99GB matched-data inputs or the MATLAB cross-spectra, which stay
on bluehive/the lab NAS per this repo's established practice).

**Final tally, 1520/1520 work units, zero errors**:

| Technique | Converged | Total | Rate |
|---|---|---|---|
| single-taper | 0 | 380 | 0.0% |
| FastMspec | 99 | 380 | 26.1% |
| Mspec | 58 | 380 | 15.3% |
| MspecBestK | 100 | 380 | 26.3% |

single-taper's 0% is a real, mechanistic result (unsmoothed coherence defeats the picker's
cycle-jump tracking entirely, not a bug) -- see the "first look at partial results" log entry.
FastMspec and single-taper are cross-validated against Sayan Swar's own precomputed MATLAB
cross-spectra where available; Mspec and MspecBestK have no MATLAB reference (never computed
upstream) and are Python-only results, flagged as such throughout.

**Each `results/<pair>__<technique>.json`** is a `WorkUnitResult` (see `work_unit.py`): distance,
convergence, best-fit template delta, the picker's own quality scalars
(`bad_quality_fraction`, `freq_coverage_fraction`, `n_candidate_crossings`, `n_accepted_picks`,
`mean_amp_ratio`), template-scan stats (`n_templates_converged`/`n_templates_scanned`), MATLAB
cross-validation where it exists, runtime, and any error. Per-template detail (all ~33 attempts,
not just the winner) and the raw coherence curve/zero-crossings were **not** captured in this
round -- see the two-round decision in the Stage 4 log for why, and Round 2's design for what a
richer capture would add.

**Used by**: the Stage 5 notebook rewrite (not yet started), and the scalar-mining analysis
already summarized in `docs/notebook5_revamp_progress.md` (resolution-bandwidth confirmation,
template-scan stability, MATLAB-mismatch-vs-convergence check).
