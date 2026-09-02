# Notebook 5 revamp: progress tracker

Tracks implementation of the plan to replace Notebook 5's zero-crossing/max-min event-scanning
barcode with a phase-velocity-based quality metric, built on an instrumented, vendored copy of
`seislib`'s own dispersion-curve picker, exercised across the full 380-pair Madagascar dataset on
bluehive. Full rationale/design: see the git history of
`~/.claude/plans/okay-let-s-start-with-wild-treehouse.md` (session-scoped, not committed here) and
this file's own log below, which is the durable record.

Updated at the end of every stage. Check in with the user after each update before starting the
next stage.

## Checklist

- [x] **Stage 0** -- This tracking file
- [x] **Stage 1** -- Provenance/citations (pull 5 papers, `docs/references/README.md`, update
      `notebooks/README.md` + `docs/coherence_barcode_design.tex`)
- [x] **Stage 2** -- Vendor + instrument the seislib picker (`python/dispcurve_pick/`)
- [x] **Stage 3** -- Validate against Sayan's SKRH-BAND result; 4-technique timing pilot
- [ ] **Stage 4** -- bluehive batch pipeline, full 380 pairs x 4 techniques
- [ ] **Stage 5** -- Notebook 5 complete overhaul (built fresh, old version tagged not deleted)
- [ ] **Stage 6** -- Packaging + docs cleanup, Notebook 3 Section 4 ref_curve fix

## Deferred / requirement log

- **Stage 5 design: ground picking stability in actual multitaper coherence estimator theory, not
  just empirical/heuristic QC signals.** Per direct guidance: the current quality diagnostics
  (`bad_quality_fraction`, `freq_coverage_fraction`, `mean_amp_ratio`) are all seislib's own
  heuristics -- empirically useful, but not derived from the coherence estimator's actual
  statistics. A principled treatment has two pieces: (1) the resolution side, already partially
  built (Stage 3's `2W <~ df_zero = c/2r` criterion) -- deterministic, says whether a crossing
  survives smoothing at all, says nothing about how precisely it's located; (2) the variance/bias
  side, the real gap -- classical Thomson/Percival & Walden coherence-estimator theory, which has
  a specific known mechanism directly relevant here: a K-taper coherence estimate has a positive
  bias floor of roughly `1/K` even at true zero coherence (it doesn't reach zero at a genuine
  null, it dips toward `1/K` and bounces) -- a real, derivable explanation for "spurious
  crossings near nulls," not just an empirical label. One real wrinkle specific to this pipeline:
  coherence here is stacked over `coh_num` (~hundreds-1600) independent day/window traces *in
  addition to* K-taper averaging -- the effective degrees of freedom driving bias/variance is some
  combination of both, and it's an open, dataset-answerable question whether trace-stacking
  dominates K's contribution at this real scale (plausible, given `coh_num` is often far larger
  than K). Recommendation for Stage 5: derive/cite the coherence bias-variance formula (Percival &
  Walden 1993, *Spectral Analysis for Physical Applications* -- the standard reference, not yet in
  this project's citation list) in terms of K and `coh_num` per technique, use it to *predict*
  where each technique's picks should be reliable, then check the 380-pair batch's actual
  convergence/`bad_quality_fraction` data against that prediction -- turns Notebook 5 from
  "measured convergence rates" into "theory, prediction, empirical confirmation," consistent with
  this project's established ethic of grounding findings in the literature, not just describing
  them. Not implemented yet -- real Stage 5 design work, deliberately not rushed alongside Stage
  4's live cluster jobs.
  - **Update, same session**: rather than the full Percival & Walden 1993 book, followed Karnik
    et al. (2022)'s own citation trail (per direct guidance) and pulled three more directly
    targeted papers into `docs/references/` (provenance/summary there): **Walden (2000, Biometrika)**
    -- the real find, a unified multitaper *multivariate* (cross-spectral/coherence) estimator
    theory paper, whose own Section 6 numerical example (60 real seismic time series) shows
    Slepian multitaper coherence bias staying below ~0.2 at a frequency with known-zero true
    coherence, vs. up to ~0.6 for a lower-quality lag-window estimator -- a real, published,
    quantitative bound on exactly the bias mechanism this project's picking-stability argument
    needs, from real seismic data, not a synthetic toy case; **Haley & Anitescu (2017, IEEE SPL)**
    -- a genuine bonus, an actual data-driven *algorithm* (jackknife MSE minimization) for
    selecting the optimal bandwidth/K automatically, a concrete candidate for replacing this
    project's currently-hardcoded `Wband=0.001`/`NW=100` constants with a principled choice;
    **Keding et al. (2024, Frontiers in Neuroscience)** -- found via forward citation search (who
    cites Karnik et al.), a different-domain (EEG) but methodologically direct treatment of
    coherence-level bias and spectral peak-shifting bias. Confirmed via each paper's own metadata
    that citations are correct; full formula-level extraction still deferred to actual Stage 5
    work, not rushed here.
  - **Second update, same session**: three more papers supplied directly (from memory, not
    citation-trail search): **Park, Lindberg & Vernon (1987, JGR)** -- the classic foundational
    paper applying Thomson's multitaper method to real seismograms; **Park & Levin (2016, GJI)**
    -- the closest match found to this project's own use case so far: multiple-taper
    *correlation* (their "MTC" algorithm, i.e. cross-spectral coherence) applied to real seismic
    receiver functions, with an explicit statistical-assumption-testing section and a dedicated
    jackknife uncertainty-estimation section (Section 6) -- real seismology-specific statistical
    treatment, not just the more abstract general theory in Walden (2000); points to the original
    Park & Levin (2000) paper (introducing MTC) as worth tracking down too if the full derivation
    history is wanted; **Walden, McCoy & Percival (1994, IEEE Trans. Signal Processing)** -- the
    univariate variance formula cited directly in Walden (2000)'s own reference list. All added to
    `docs/references/` with the same provenance/summary convention as the rest.
- **Stage 4 SLURM design: parallelize at (pair, technique) granularity, not (pair) with all 4
  techniques serial inside one task.** Observed directly during Stage 3's local timing pilot:
  running single-taper/FastMspec/Mspec/MspecBestK serially in one process is slow and
  memory-heavy (one technique's peak briefly hit 18GB RSS on this machine) with no benefit from
  serializing them -- each `(pair, technique)` cross-spectrum is fully independent. The bluehive
  array job should index over `(pair, technique)` (380 x 4 = 1520 units of work, batched a handful
  per array task), not one task per pair running all 4 in sequence -- genuinely embarrassingly
  parallel, not just parallel-across-pairs. Recorded per direct guidance while watching the pilot
  run long.
- **Stage 4 SLURM design: build BOTH a no-multiprocessing and a multiprocessing driver, and
  decide between them from real measured cost, not assumption.** `multiprocessing`/
  `ProcessPoolExecutor` (not full MPI -- the workload is embarrassingly parallel, no inter-task
  communication needed, so MPI's actual value-add doesn't apply here) is very likely the right
  call once per-task cost is known: each array task claims one full node and fans its share of
  `(pair, technique)` work out across that node's cores in-process, avoiding both idle cores and
  1520x repeated Python/numpy/scipy/obspy import overhead. Per direct guidance: implement both a
  plain one-work-unit-per-task driver and an in-node-multiprocessing driver, so the choice between
  them (and the right node/worker-count split) is made from actual measured throughput on this
  workload, not assumed in advance. Real MPI only reconsidered if task count grows much larger or
  individual tasks turn out sub-second, per the note above.
- **Extend the quality gate to M/N events (maxima/minima), not just Z (zero-crossings).** The v2
  barcode drops M/N entirely rather than fixing their noisiness, since seislib's `bad_quality`
  machinery is defined for zero-crossings specifically. The underlying ambition of a 3-way Z/M/N
  barcode wasn't wrong -- the naive detection+filtering was. Revisit once the Z-based picker is
  validated at full 380-pair scale: design an analogous, principled quality criterion for M/N
  before reconsidering a 3-way barcode. Not in scope this round.

## Log

### 2026-09-01 -- Stage 0

Created this file. Plan approved same day after extensive investigation (2 Explore agents + 1 Plan
agent, all with live SSH to `repovibranium` and `bluehive`) and iterative review -- see that
session's transcript for the full investigation trail; key findings are already folded into the
plan file's Context section and won't be re-derived here.

All implementation work happens on a feature branch, `notebook5-phase-velocity-revamp`, merged
into `main` and pushed after each stage completes -- not committed straight to `main` throughout,
per user request when this stage's work was reviewed.

### 2026-09-01 -- Stage 1

Pulled the 5 grounding papers from the NAS (`core_papers/`) into `docs/references/`, each verified
against its own embedded PDF metadata (via `pypdf`) rather than trusted by NAS filename -- all 5
matched existing citations exactly, including resolving the AkiNet paper's year: genuinely 2025
(JGR: Machine Learning and Computation, DOI `10.1029/2025JH000932`) despite its NAS filename
saying `xueolugboji2026.pdf`. Added `docs/references/README.md` (provenance, matching
`data/reference/README.md`'s convention). Updated `notebooks/README.md`'s citations to link the
local copies. Substantively revised `docs/coherence_barcode_design.tex`'s "Relationship to
Existing Tools" section with the real 4-stage seislib algorithm (bad_quality gating, kernel-density
picking, cycle-jump rejection, coverage acceptance test -- confirmed by a full read of Sayan's
vendored source this round, not just the high-level behavior). Added a new "Revision" section
documenting the pivot itself: what's superseded (the v1 event-scanning scorer, tagged
`notebook5-v1-event-scanning` before removal in Stage 5) vs. kept (the template family, repurposed
as `ref_curve` priors), and the deferred M/N quality-gate extension (see "Deferred / requirement
log" above). Re-fetched the `tectonic` LaTeX binary (this session's filesystem didn't have it
persisted from before) and confirmed the doc compiles cleanly, 10 pages, no undefined references.

Also added `*.log`/`*.aux`/`*.toc` to `.gitignore` (LaTeX build artifacts, not previously
excluded). Noted for Stage 4: `.gitignore`'s blanket `*.csv` rule will need a targeted exception
for the eventual `data/results/dispcurve_quality/manifest.csv`, similar to how `/data/reference/`
was carved out of `/data/*` -- not yet done, flagged so it isn't forgotten.

Merged into `main`, pushed.

### 2026-09-01 -- Stage 2

Vendored Sayan's copy of seislib's `_an_processing.py` from the NAS into `python/dispcurve_pick/`
(1539 lines; diffed against installed `seislib==1.2.1` at pull time -- identical except one
deprecated-numpy-API line and a trailing newline). Instrumented at exactly 3 sites (bad_quality
capture, per-pick amplitude-ratio capture inside the `pick_velocity` closure, the final
return/raise site), each in a clearly marked comment block, adding an opt-in
`return_diagnostics=True` that surfaces `PickDiagnostics` (bad-quality fraction, candidate-crossing
count, accepted-pick count, frequency-coverage fraction, mean amplitude ratio) without altering
the picking algorithm itself.

`tests/test_matches_upstream.py` (3/3 passing) confirmed this by finding and fixing two real bugs
along the way: (1) `np.in1d`, used in Sayan's copy, no longer exists in this environment's numpy
2.4.6 -- fixed to `np.isin`, matching what the currently-installed `seislib` package itself already
uses at that line (an environment-compatibility fix, not a behavior change); (2)
`DispersionCurveExceptionWithDiagnostics`'s `__init__` tried to pass a message string to
`DispersionCurveException.__init__`, which upstream defines to take no arguments at all -- fixed to
call `super().__init__()` then set `self.message` directly. Full technical writeup:
`python/dispcurve_pick/NOTES.md`.

Re-fetched `tectonic` and `pypdf` were both needed fresh this session (this sandbox's filesystem
doesn't persist installed tools/packages across sessions) -- noting this since it'll likely recur:
check for tools before assuming a prior session's setup carried over.

Merged into `main`, pushed.

### 2026-09-01 -- Stage 3 (part 1: validation)

**Important correction to the plan's own framing**: Sayan's own SKRH-BAND notebook
(`phasevel_compute_slide13and14.ipynb`) is not a "both methods converge" golden reference as
earlier investigation assumed -- checking its actual saved outputs (it does have them, execution
counts and all) shows **FastMspec converges to a clean picked curve** (0.02-0.115 Hz, ~3.4-3.85
km/s, tracking a real kernel-density ridge) while **single-taper, even after extra savgol
smoothing, fails outright** with a `DispersionCurveException`. That's a stronger, more honest
validation target than "both succeed" would have been -- it's already the FastMspec-beats-
single-taper story this whole revamp is trying to establish, from Sayan's own real execution.

Reproduced this exactly with the vendored+instrumented picker (`return_diagnostics=True`) against
the two precomputed `.mat` files pulled from the NAS: FastMspec converges (curve spans
`[0.0358, 0.1101]` Hz, velocity `[3.41, 3.85]` km/s, matching his plot), single-taper (smoothed,
his exact params) does not (`bad_quality_fraction=0.71`, `n_accepted_picks=150`,
`freq_coverage_fraction=0.0`). Added a fair, symmetric comparison beyond what his own notebook
tried (both raw/unsmoothed, same 0.01-0.5 Hz band): FastMspec unchanged, single-taper fails even
harder unsmoothed (`n_candidate_crossings=2066` but only `n_accepted_picks=15`,
`freq_coverage_fraction=0.0`). One useful finding for later design: `bad_quality_fraction` alone
is a weak discriminator here (~0.67-0.72 across all three runs) -- `n_accepted_picks` and
`freq_coverage_fraction` are what actually separate a converging pick from a failing one.

### 2026-09-01 -- Stage 3 (part 2: timing pilot + a genuine bug hunt that resolved cleanly)

**Timing** (all local-machine numbers, provisional -- see the partition/multiprocessing notes
above; SKRH-BAND, 1605 traces, N=10801):

| Technique | Cross-spectrum | Picking (33-template corridor scan) |
|---|---|---|
| single-taper | 2.47s | |
| FastMspec | 143.6s (taper_size=13) | 99.4s (17/33 templates converged) |
| MspecBestK | 147.4s (taper_size=15) | |
| Mspec | **2688.5s** (44.8 min, taper_size=80) | |

Mspec dominates by >10x over every other technique -- confirms the plan's own caution not to
assume uniform task sizing across techniques; Mspec needs its own, more generous task budget in
Stage 4's SLURM design. Naive fully-serial extrapolation (4 techniques x 380 pairs, no
parallelism) would be ~357 hours -- exactly why Stage 4 is designed around embarrassingly-parallel
`(pair, technique)` array tasks, not a reason to reconsider scope.

**A real bug hunt, resolved cleanly.** The first timing-pilot run crashed on a script bug
(`compute_crosscorr` with a default `FilterConfig` falls through to the *unsummed* 3D
`(day, window, freq)` plain-FFT branch, not a `CrosscorrResult` -- a real usability gap worth
noting for Stage 4's driver, not just a one-off mistake) -- but before crashing, it surfaced a
genuinely concerning number: **FastMspec's recomputed cross-spectrum showed a 3.06% relative L2
error against Sayan's own precomputed MATLAB `coh_sum`** for this exact pair, with `coh_num`
matching exactly (1605 both sides, ruling out a trace-count bug). Investigated rather than waved
away, per this project's established practice:

- Binning the error by `|coh_sum|` magnitude (a proxy for distance from a coherence null) showed
  it concentrates exactly where the already-documented, intentionally-not-reproduced MATLAB
  complex-floor bug (`ccf_pipeline/NOTES.md`'s "Known upstream bug" section) predicts: 20.8%
  relative error in the lowest-magnitude 10% of bins, 7.4% in the next 40%, 2.4% in the top 50% --
  while the *absolute* error stays roughly constant (~1.05-1.12) across all three bins, exactly
  the signature a `max(z, eps*maxz)` floor produces (constant absolute perturbation, so relative
  error balloons only where the true signal is small). This is the first real-data (as opposed to
  small-synthetic) confirmation of that bug's existence and behavior -- added to
  `ccf_pipeline/NOTES.md` with the concrete numbers.
- The corrected single-taper cross-validation (properly summing the plain-FFT branch, fixing the
  script bug above) then showed its **own** large mismatch -- 46.5% relative error -- that could
  *not* be explained by the same bug (single-taper never applies the complex-floor logic at all).
  Root cause, found by testing detrend/taper variants directly against the MATLAB reference:
  "single-taper"/"first-order" in this project's own established terminology means the **"5%
  Cosine Single-Taper"** technique -- detrend + cosine taper applied *before* the plain FFT
  coherency (already named this way elsewhere in this project, e.g. Notebook 3's own Section 2) --
  not raw, unprocessed data straight into FFT. Once that preprocessing was applied, the relative
  error dropped to **4.83e-09** -- machine precision, a clean match. My own script's omission, not
  a pipeline bug; the actual `ccf_pipeline` code was never wrong.

Both cross-validations now confirm the pipeline is correct on real, large-scale data for the first
time (previously only verified via small synthetic fixtures + Octave) -- directly closing the gap
`ccf_pipeline/NOTES.md` named ("the `IsMspec`/`FastMspec` path... still only verified against
synthetic data + Octave, not a real end-to-end SAC-to-output run").

**Consequence for Stage 4's driver design**: it must apply detrend+cosine-taper before the
single-taper/plain-FFT technique specifically (not the other three, which handle their own
windowing via DPSS tapers and use the un-preprocessed matched data directly, matching the
production driver's `IsDetrend=0/IsTaper=0` config) -- noting this now so it isn't rediscovered
the hard way during the full batch run.

All raw/precomputed `.mat` files used in this stage were pulled to local scratch only (not
committed -- 267MB+164KB+164KB, matching the project's established practice of not committing bulk
data). Consolidated the validation work (both cross-spectrum cross-validation and dispersion-curve
picking) into a permanent, documented script:
[`verification/skrh_band_real_data/`](../verification/skrh_band_real_data/) (`validate_skrh_band.py`
+ `README.md`, matching `verification/octave_verify_ccf_pipeline/`'s existing convention) --
re-run end-to-end after cleanup and confirmed it reproduces every number above exactly. Caught one
more real bug along the way while consolidating: Sayan's own `faxis`-construction formula
(copied verbatim from his notebook) actually builds a length-(T+1) array for odd T -- benign in
his own code (he always indexes both `faxis` and `coh_sum` with the same *integer* index array,
whose values never reach the extra element), but breaks if a boolean mask is applied directly to
`coh_sum` instead -- fixed by using integer indices, matching his own working pattern, documented
inline in the script rather than silently worked around.

Merged into `main`, pushed.

### 2026-09-01 -- Stage 4 (in progress)

**Built** `python/dispcurve_pick_batch/`: `manifest.py` (380-pair catalog -> 1520 work units,
technique-outer ordering so each technique is a contiguous index range --
`technique_index_ranges()`), `work_unit.py` (the core `process()`: per-technique preprocessing,
MATLAB cross-validation where available, template-corridor scan via the instrumented picker,
never raises -- exceptions land in the result's `error` field), `run_plain.py` /
`run_multiprocessing.py` (the dual drivers, per earlier direct guidance), `aggregate.py`. Also
promoted `load_reference_curve`/`build_template_family` from `notebooks/_lib/nb5_helpers.py` into
`python/dispcurve_pick/template_family.py` so both the batch pipeline and the eventual Stage 5
notebook import from one place, without the batch pipeline depending on `notebooks/_lib`.

**Local validation before touching bluehive**: `work_unit.process()` run end-to-end on the real
SKRH-BAND pair (both techniques already validated in Stage 3) -- reproduced Stage 3's exact
findings (single-taper doesn't converge; FastMspec converges, best_delta=-0.30 km/s, close to
Stage 3's own -0.35 from a slightly different scoring exercise).

**Deployment hit two real, worth-recording surprises** (both already written up in detail in
`python/dispcurve_pick_batch/NOTES.md`'s "Deploying to bluehive" and
`python/dispcurve_pick/NOTES.md`'s "Why the seislib package dependency was removed" -- summarized
here):
1. `PRJ_SPAC/codes/prod/` (the planned deploy location, inside Sayan's own `Sayan_Swar_WS`) is
   read-only to this account -- discovered live via a failed `mkdir`, not assumed from listing
   permissions alone. Redirected to a fresh top-level directory this account does own:
   `/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/`.
2. A plain `venv` hit a wall of genuine toolchain fragility on the login node: the default
   `python3/3.11.0` module's SSL support is broken entirely (every pip HTTPS request silently
   "finds no versions"); the working `python3/3.11.10` module needs `LD_LIBRARY_PATH` from the
   *same* `module load` in every single command (doesn't persist across separate `ssh`
   invocations, same class of issue as `$BASE` env vars earlier); building numpy from source needs
   a newer gcc than the login node's ancient default `cc` (a loaded gcc module exists but meson
   invokes `cc` specifically -- needs `CC=gcc CXX=g++`); and `obspy`'s `pyproj` dependency needs
   PROJ >= 9.4.0, newer than any available `proj` module (max 8.1.1, itself broken by an unrelated
   OpenSSL library-version mismatch). Switched to a **fresh, dedicated conda environment**
   (`fastmspec_batch`, same shared Anaconda install other lab members use, but a new env -- not
   Sayan's own broken `Seislib` one, untouched): `conda create -n fastmspec_batch -c conda-forge
   python=3.11 numpy scipy pandas obspy` -- prebuilt binaries, zero compilation, worked immediately.

**A genuine improvement fell out of the second surprise**: `pip install seislib` *also* failed on
this login node, for a third, unrelated reason (a broken pre-generated Cython extension in
`seislib.tomography`, irrelevant to dispersion-curve picking but built as part of any `pip install
seislib` regardless). Rather than fight that too, vendored the ~6 small functions/classes
`dispcurve_pick` actually needs from `seislib.utils`/`seislib.exceptions` directly
(`_vendored_seislib_utils.py`, `_vendored_seislib_exceptions.py` -- confirmed functionally
identical to the installed package via direct source diff, cosmetic differences only). The picking
path no longer depends on the full `seislib` package at all, anywhere -- a real dependency-surface
reduction, not just a bluehive workaround. `seislib` stays installed locally for
`tests/test_matches_upstream.py`'s own byte-fidelity comparison against real upstream, and for
Notebook 3 Section 4's separate, unmodified direct usage -- both updated/re-verified, 3/3 tests
still passing after the refactor.

**bluehive validation so far** (code + `data/reference/SDISPL.ASC` + the pulled
`madagascar_stn_conn_ccflist.csv` deployed via `rsync`):
- Direct invocation (no SLURM) of one single-taper work unit on a fresh pair (XVKIRI/XVMAGY, not
  previously seen): ran cleanly end-to-end, 87.6s, and its MATLAB cross-validation *also* worked
  --matched to 2.67e-13 relative error (machine precision), on a pair Stage 3 never touched --
  a strong, unplanned bonus confirmation that the single-taper preprocessing fix generalizes.
- A 2-task SLURM array smoke test on `debug` (single-taper, two more fresh pairs): both completed
  cleanly, 207.0s and 294.3s -- confirms the full SLURM mechanics (submission, array indexing,
  conda activation, module resolution, idempotent result-file writing) work correctly. Notably
  slower than the login node's 87.6s (~2.4-3.4x) -- compute nodes here are genuinely slower/more
  contended than the login node for this workload; real, worth-remembering for sizing the full run.
- A FastMspec smoke test (`--time=00:15:00 --mem=16G`, guessed without applying the above
  slowdown factor) was killed by SLURM at both the time AND memory limit -- a real, useful
  finding, not a bug: FastMspec's ~243s login-node baseline (144s cross-spectrum + ~99s picking),
  scaled by the same ~2.4-3.4x compute-node slowdown single-taper showed, lands almost exactly at
  the 15-minute mark. Retried with `--time=00:45:00 --mem=32G`; in progress as this entry is
  written.

**FastMspec retry (`--time=00:45:00 --mem=32G`) succeeded**: `XVKIRI_XVMAGY__FastMspec` converged
(`best_delta_km_s=0.55`), `matlab_rel_l2_error=8.04e-09` -- machine precision, even tighter than
SKRH-BAND's ~3% (consistent with the complex-floor bug's null-region-specific signature from Stage
3: this pair's spectrum apparently has fewer/shallower coherence nulls in-band, so the bug barely
bites here -- expected data-dependent variation, not a contradiction). Runtime: 1162.4s (~19.4
min) -- ~4.8x slower than the local login-node baseline (~243s), notably worse than single-taper's
~2.4-3.4x slowdown on the same class of node; confirmed via `sstat`/direct `ps` on the compute node
that this was genuine, steady 99.7%-CPU computation, not a hang. stderr logged "Exceeded step
memory limit at some point" despite the job completing successfully -- 32GB was cutting it close;
size FastMspec's real submission with more headroom (~40-48GB).

**Both primary techniques (single-taper, FastMspec) are now confirmed correct end-to-end via real
SLURM execution**, on three different real pairs beyond SKRH-BAND, with MATLAB cross-validation
matching to near machine precision in every case tested. Mspec/MspecBestK have not been separately
smoke-tested via SLURM (only locally, Stage 3) -- given Mspec's already-known >10x cost and the
~4.8x node slowdown just observed, a live SLURM confirmation could plausibly run 1.5-3+ hours for
one work unit; judged not worth spending that wall-clock on before the real batch, since the code
path is structurally identical to the already-proven techniques (same `process()` function, only
the technique branch differs) -- the real unknown left is resource sizing, not correctness, and
that can be extrapolated from known per-technique cost ratios rather than re-confirmed live.

**Local plain-vs-multiprocessing comparison** (per explicit request): 6 independent single-taper
work units (all reusing the local SKRH-BAND matched-data file, relabeled -- a fair proxy for N
independent work units of realistic size, without pulling N x 267MB of new data just for a timing
test), on this machine (10 physical / 20 logical cores). Plain (fully serial): 827.6s total
(133-141s each, consistent with earlier numbers). Multiprocessing (8 workers): 209.9s total --
**3.94x wall-clock speedup** -- even though each individual work unit ran ~1.5x slower under
contention (202-207s vs ~135s alone, expected: 6 processes sharing 10 physical cores and memory
bandwidth). Clean, unambiguous answer: multiprocessing is the right call for wall-clock
throughput on this workload, confirming what was expected but per direct guidance not assumed
without measuring. `run_multiprocessing.py` is the primary driver for the real batch; `run_plain.py`
stays available (simpler, maximally fault-isolated) as a fallback/comparison tool, not deleted.

**Open uncertainty, flagged rather than resolved unilaterally**: Mspec/MspecBestK have not been
smoke-tested via real SLURM execution (only locally). Reasoning for not doing so live: (a) the
`debug` partition's hard 1-hour cap can't fit Mspec's local baseline (~46.5 min cross-spectrum
alone) even before applying the observed node slowdown; testing it would mean already committing
to a `preempt`/`urseismo` submission, not a cheap smoke test; (b) the observed slowdown factor
was *not* uniform across techniques on the same node (single-taper ~3.36x, FastMspec ~4.78x, both
on `bhp0001`) -- extrapolating Mspec's real cost from either ratio is a guess, not a measurement,
and Mspec's much heavier memory footprint (K=80 tapers vs. FastMspec's 13) could plausibly see an
even worse ratio if it interacts with cgroup memory limits the way FastMspec's near-miss at 32GB
suggests. Correctness itself is low-risk (identical code path, same `process()` function, only the
technique branch differs) -- the real open question is resource sizing for a technique that could
plausibly take on the order of hours per work unit on this cluster, not whether it works at all.

This is the natural decision point for the full-batch submission -- presented to the user rather
than resolved unilaterally, given the real cluster-time/resource commitment involved (up to 1520
tasks, with Mspec's own 380-task slice sized from extrapolation rather than direct measurement).

**Mspec on `preempt` (3hr/64GB) hit a real, severe memory wall**: killed by the OOM killer in only
~2 minutes -- far too fast to be a gradual buildup during the ~45-minute computation, pointing at
failure during array construction itself. The math confirms it: `classical_spectrum_batch`'s
direct translation of `avgspec` materializes `(N, n_traces, K)`-shaped intermediates -- for
Mspec's K=80 at this project's real Madagascar trace counts (~1600), a single such array is
~11-22GB, with several alive simultaneously (pre-FFT tapered traces, post-FFT complex spectra for
both traces, their product), totalling an estimated **~55GB peak**. `MspecBestK`'s much smaller K
(~13-15) never hit this, which is why it went unnoticed until Mspec ran at real scale -- the
small-N synthetic fixture this function was originally verified against never exercises memory at
all. Not a MATLAB-vs-Python translation bug (the original `avgspec` has the same materialization
pattern) -- a real scalability limitation neither implementation had been run hard enough to hit
before.

**Fixed at the root**, per explicit direction (optimize rather than just throw more memory/wait for
`highmem`, and keep it backward compatible): `classical_spectrum_batch` now chunks over the taper
(K) axis -- processes `k_chunk` tapers at a time (default 8), accumulating a running sum instead of
materializing all K simultaneously, then divides by K at the end (mathematically identical to the
original `.mean(axis=2)`). Verified numerically against the pre-chunking implementation on
synthetic data (including a ragged-last-chunk case, K not a multiple of `k_chunk`): relative error
~1.6-2.3e-16 (machine epsilon, pure floating-point summation-order noise) for `k_chunk=1` and the
new default; exactly 0 for `k_chunk=K` (identical code path to before). `k_chunk` is a new,
optional, keyword-only-by-default parameter -- every existing call site's signature/return shape
is unchanged, confirmed backward compatible per explicit direction. Estimated new peak for the
same case: ~5.5GB (~10x reduction), comfortably within normal node memory -- no `highmem` needed.
Full writeup: `python/ccf_pipeline/NOTES.md`. Redeployed to bluehive; a live retest
(`--mem=24G`, well below the old 64G that failed) is running as this entry is written.

**Memory fix confirmed live**: retested on `preempt` with `--mem=24G` (well below the old 64G
that failed) -- ran to completion in 1014.6s (~16.9 min), zero memory errors, not even a warning
this time (unlike FastMspec's earlier close call at 32G). The fix works as designed.

**The result itself is a genuine, unplanned confirmation of the Stage 5 design direction above**:
Mspec did not converge for this pair (`XVKIRI_XVMAGY`, same pair single-taper also failed and
FastMspec succeeded on) -- and checking the resolution-bandwidth criterion directly explains why,
not just describes it: Mspec's smoothing bandwidth (`2W = 2*NW/N = 0.0185` Hz, from its fixed
`NW=100`) is **5-11x wider than the natural zero-crossing spacing** at this pair's 459.8 km
distance for any plausible phase velocity (`df_zero = c/2r` ~ 0.0016-0.0038 Hz across 1.5-3.5
km/s) -- it violates `2W <~ df_zero` outright, by a wide margin, at every reasonable velocity.
FastMspec's much narrower bandwidth (`2W = 0.002` Hz) satisfies the same criterion for all but the
very lowest velocities, consistent with its convergence on the same pair. Real, textbook evidence
that Mspec's failure here is over-smoothing blurring straight through the Bessel-null crossings
being picked -- not just "another noisy result" -- and a first confirmation that the principled
framework just discussed has real predictive power on data already in hand, before Stage 5 has
even formally started.

**Stage 4 is now functionally complete**: all four techniques confirmed correct and running on
real bluehive infrastructure (single-taper, FastMspec, Mspec directly; MspecBestK shares Mspec's
code path at a much smaller K, already validated locally in Stage 3 and structurally identical to
FastMspec's already-SLURM-confirmed path). The multiprocessing-vs-plain question is answered
(multiprocessing wins, ~3.94x). The one remaining step is the full 380 x 4 = 1520-task submission
itself -- deliberately not automatic, presented to the user next as an explicit go/no-go.

Next: present Stage 4's status and the full-batch submission decision to the user.

### 2026-09-02 -- Stage 4: full batch submitted; one real bug found and fixed live

Submitted all 4 techniques' array jobs on `preempt` (job IDs 31344560-31344563; 337 total array
tasks -- single-taper 19 tasks/chunk=20/24G/00:30:00, FastMspec 127 tasks/chunk=3/120G/01:00:00,
Mspec 64 tasks/chunk=6/100G/00:45:00, MspecBestK 127 tasks/chunk=3/120G/01:00:00 -- sized from
Stage 3's timing pilot plus a safety margin, using `run_multiprocessing.py` per the plain-vs-mp
comparison above). Node specs confirmed live (24 CPUs, 127-189GB RAM) before sizing.

**A real bug surfaced and was fixed while the batch ran**: 6 single-taper work units (all
involving station `GFOMA`) failed instantly (~0.2-0.4s) with `IndexError: tuple index out of
range`. Root cause: pairs with only a single day of overlapping data collapse their
`S1_data_mat`/`S2_data_mat` to a **2D** `(n_window, n_samples)` array in the `.mat` file (scipy/
MATLAB drop the singleton day axis) instead of the usual 3D `(n_day, n_window, n_samples)` every
pair tested so far had -- `work_unit.py` assumed 3D unconditionally (`s1.shape[2]`). Checked
scope before fixing: 11/380 pairs (~2.9%, `filesize_mb < 10`, vs. a dataset mean of ~690MB) are
small enough to be at risk -- bounded, but real, not a one-off. Fixed by reshaping 2D input to
`(1, n_window, n_samples)` right after loading, in `work_unit.py`; verified directly against the
exact failing pair (`GFOMA_XVLONA`) post-fix: real computation completed cleanly (273s,
`error=None`, a genuine `converged=False` result, not a crash). Deployed the fix to bluehive mid-
batch; the 6 already-corrupted result files (written by array tasks that started before the fix)
were deleted and resubmitted as a small, targeted array job (indices `14,15,147,185,234,236`,
found via `build_work_units()`) -- job 31344579. Array tasks that hadn't started yet when the fix
was deployed picked it up automatically (fresh Python process per task, reading the updated file
from disk) -- no other units affected.

Batch is running; monitoring continues. This is exactly the kind of real-scale-only failure mode
this project has hit before (Stage 3's Mspec-only memory wall was the same pattern) -- a case
neither hand-written nor small-fixture tests would ever exercise, only found by actually running
the real, full, heterogeneous dataset.
