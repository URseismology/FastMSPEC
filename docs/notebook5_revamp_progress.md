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

- **Stage 5 design, sharpened framing (2026-09-02, for next session): `K`/`W` needs an upper bound
  AND a lower bound, and only the upper bound is principled so far.** The resolution-bandwidth
  criterion (`2W <~ Delta_f_zero = c/2r`, below) gives the **upper bound**: `W` (hence `K`, via
  `K ~ 2NW-1`) too large over-smooths and merges/misses real zero-crossings -- this is what's
  already confirmed against real Round-1 convergence-vs-distance data. Still missing: a **principled
  lower bound** on `K` -- how small can `K` go before the coherence estimator's own variance/bias
  (the `~1/K` bias floor near true-null coherence, Walden 2000) starts *manufacturing* spurious
  crossings instead of merely missing real ones. Today's `K` values (13-15 MspecBestK, 80 Mspec,
  FastMspec's fused equivalent) are hand-picked, not derived from either bound. This is the first
  thing to review next session, once Round 1's results are in hand to check any candidate lower-
  bound criterion against real convergence/bad-quality-fraction data, before Round 2 is designed.
  - **Update (2026-09-03): the upper/lower-bound derivation and the memory-architecture argument
    are coupled, not independent -- and this decouples Stage 5's bandwidth search from any memory
    concern before it's even done.** Per direct guidance: whatever principled `W` eventually falls
    out of the resolution upper bound and the (not yet derived) bias-variance lower bound, it will
    land somewhere between them -- and separately verified (`FastMultitaper` run directly across
    NW=5 to NW=400 at N=10801) that FastMspec's transition-region correction count `r` stays
    essentially bounded (11-20) across that *entire* range while `K` grows linearly the whole way
    (5 to 790) -- see the K-chunking writeup below for the full table. So this isn't "FastMspec
    happens to be cheap at today's ad-hoc `Wband=0.001`" -- it's "FastMspec will be cheap at
    whatever `W` Stage 5's principled search lands on," because `r`'s boundedness is a DPSS
    transition-region property independent of where in that range `NW` falls. Practical
    consequence for Stage 5: the bandwidth search can be run purely on statistical grounds
    (resolution vs. bias-variance), without needing to separately re-validate memory/chunking
    strategy for whatever `W` it recommends -- FastMspec absorbs that automatically, which
    classical Mspec/MspecBestK structurally cannot (their memory scales with `K` directly, so a
    revised `W` mid-design would mean a revised memory budget too). This is now the headline
    architectural argument for FastMspec in Stage 5, not merely a secondary efficiency note.
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

The GFOMA fixup (job 31344579) is now confirmed complete: all 6 targeted array tasks (indices
14, 15, 147, 185, 234, 236) show `COMPLETED, 0:0` in `sacct`, runtimes 4m22s-20m40s, all with real
`error=None` results -- the 2D-shape bug is fully resolved across the dataset.

### 2026-09-02 -- Stage 4: a second real bug -- FastMspec per-pair runtime variance blows the
1-hour task time limit

While monitoring the batch, found `XVANTS_XVLONA__FastMspec` (array index 388) and
`XVBITY_XVLAHA__FastMspec` (index 405) missing after their tasks (job 31344561, array indices 2
and 8) hit `TIME LIMIT`/`CANCELLED` in their `.err` logs. Checked which of each chunk-of-3's units
had already written results (via `build_work_units()` to map array index -> work-unit indices,
then `ls results/`): 4 of the 6 covered units existed, so only these 2 specific units were
actually lost, not the full 6.

Tried `scontrol update JobId=31344561 TimeLimit=02:00:00` (and the same for 31344563/MspecBestK,
same risk) to protect everything still in flight before it hit the same wall. This worked for
pending-but-not-yet-started array indices (now at a 2-hour cap), but returned "Access/permission
denied" for the 6 tasks already running at that moment (indices 12, 13, 15, 17, 18, 20) --
apparently SLURM here won't let a user raise the time limit on an array task once it's actively
running, only while still pending. Those 6 stayed capped at the original 1 hour and were flagged
as still at risk. Submitted a first fixup for the 2 confirmed-lost units immediately
(`sbatch --partition=preempt --array=388,405 --time=02:00:00 --mem=120G --cpus-per-task=1
submit_plain.sbatch` -> job 31344667).

**The 6 flagged-at-risk tasks did then time out**, confirmed via `sacct -j 31344561`:

```
31344561_12     TIMEOUT      1:0   01:00:29
31344561_13     TIMEOUT      1:0   01:00:29
31344561_15     TIMEOUT      1:0   01:00:29
31344561_17     TIMEOUT      1:0   01:00:29
31344561_18     TIMEOUT      1:0   01:00:29
31344561_19   COMPLETED      0:0   00:54:33
31344561_20     TIMEOUT      1:0   01:00:01
```

(index 19 finished cleanly at 54m33s -- under the cap by a real but thin margin.) Checked each of
the 6 timed-out tasks' chunk-of-3 work units against `results/`: this time **all 18 were lost**,
none of the three per task had finished before the kill (`grep -c 'converged=' logs/mp_31344561_*.
out` returned 0 for all 6) -- worse than the first pair, not better.

**Root cause, now well evidenced rather than guessed**: FastMspec per-pair runtime has real,
large variance across the dataset, not the roughly-uniform cost the original chunk=3/1-hour
sizing assumed. Task 19's own log (the one that *did* finish) shows the spread directly:

```
XVMAGY_XVZAKA__FastMspec: converged=True  runtime=996.5s
XVLAHA_XVVATO__FastMspec: converged=False runtime=1040.5s
GFOMA_XVKIRI__FastMspec:  converged=True  runtime=3269.9s
```

One pair alone took 3269.9s (~54.5 minutes) -- close to the entire 1-hour task budget by itself,
running as just one of 3 concurrent workers in that task. Any chunk unlucky enough to contain even
one pair in this range, especially under any node contention on a preemptible partition, has
essentially no margin left for its other 2 units. Sizing chunk=3/01:00:00 from Stage 3's pilot
(which only sampled a small, likely non-representative slice of pair distances/sample counts) was
too optimistic for the tail of this distribution -- a genuine sizing lesson, not a code bug this
time.

**Fix applied**: resubmitted all 18 lost units individually (one work unit per task, `submit_plain.
sbatch`, not the multiprocessing driver, so no other unit's slowness can cost a fast unit its own
budget) with a much larger margin: `--time=02:00:00 --mem=120G --cpus-per-task=1`, array indices
`416,417,418,419,420,421,425,426,427,431,432,433,434,435,436,440,441,442` -> job 31344723.
Combined with the earlier job 31344667 (388, 405, same generous settings), all 20 known-lost
FastMspec units now have a fixup in flight. Both jobs were still `PENDING` (queued behind
`preempt`'s `JobArrayTaskLimit`/priority throttling) as of this log entry -- not yet confirmed
complete.

**Follow-up implication for Stage 5, noted but not acted on yet**: this per-pair runtime spread is
itself a data point about the dataset (larger/more-windowed pairs cost much more to compute), not
just an infrastructure nuisance -- worth a passing mention in the notebook's honest-discussion
section rather than being silently absorbed into "the batch eventually finished."

Both FastMspec fixups (job 31344667: indices 388, 405; job 31344723: the 18 indices from the
6-task loss) are now confirmed `COMPLETED` -- all 20 previously-lost FastMspec units recovered.

### 2026-09-02 -- Stage 4: a third bug -- Mspec's *entire* technique was mis-sized, not just its
tail (100% task failure), plus a real SLURM array-index limit discovered along the way

While the FastMspec fixups above were in flight, checked on Mspec (job 31344562, the classical
K=80 technique) and found something categorically worse than the FastMspec case: **every single
attempted task had timed out.** `sacct -j 31344562` showed indices 0-28 (29 tasks) all `TIMEOUT`
at almost exactly 45:xx minutes -- the task's own `--time=00:45:00` limit -- and the two still
`RUNNING` (29, 30) went on to time out identically a few minutes later. Checked `results/` against
every one of the 31 attempted tasks' work-unit ranges (`[760,766)` through `[934,946)`, from each
task's own `Processing N work units [...)` log line): **zero results existed.** Not one of the 31
attempted tasks' 186 work units had finished, let alone the untried remainder -- effectively 0/380
Mspec pairs done.

**This is not tail variance like the FastMspec case -- it's the whole technique undersized.** Stage
3's own timing pilot had already measured this, in fact: `2688.5s` (44.8 minutes) for a *single*
Mspec cross-spectrum computation on one pair (SKRH-BAND), recorded plainly in this file's own
Stage-3 log table. A task budgeted for `00:45:00` running *6 such pairs concurrently* (chunk=6,
`run_multiprocessing.py`) was never going to finish even its fastest unit with real margin, let
alone all 6 under node contention -- that pilot number should have driven the original sizing and
didn't. (For comparison, the plain-driver script's own header comment already had a more
realistic per-pair Mspec suggestion, `--time=02:00:00`, sitting right there unused -- the actual
submission for Mspec used the multiprocessing driver with a much tighter, inconsistent budget.)

**Fix, part 1 -- stop the bleeding.** Cancelled the 33 still-`PENDING` Mspec array indices
(`scancel 31344562_[31-63]`) before they repeated the identical failure; let the 2 already-running
tasks finish naturally (they too then timed out, as expected).

**Fix, part 2 -- resubmit properly: individual work units (plain driver, no shared-node
contention), a real time margin, via `submit_plain.sbatch`.** Attempting this surfaced a fourth,
independent, previously-unknown constraint: `sbatch --array=760-1139 ...` failed outright with
`"Invalid job array specification"` -- turned out to be this cluster's `MaxArraySize=1001`
(`scontrol show config`), which caps the raw array *index value*, not the count of tasks. Mspec's
real global work-unit indices (760-1139) and MspecBestK's (1140-1519) both exceed that on their
own -- a constraint the original plan never anticipated (it only ever indexed small technique-
local ranges via `run_multiprocessing.py`'s own internal offsetting, so this never surfaced until
the plain driver was pointed at Mspec/MspecBestK's real global range directly).

Fixed generally, not just worked around for this one case: added an `IDX_OFFSET` environment
variable to `submit_plain.sbatch` (default `0`, backward compatible -- every prior single-taper/
FastMspec submission this session used raw indices under 1001 and is unaffected), added to
`$SLURM_ARRAY_TASK_ID` before being passed to `run_plain.py`. Also corrected the script's own
header-comment usage example, which had been silently wrong since it was written (claiming a
`submit_plain.sbatch Mspec` positional technique argument that `run_plain.py` never accepted --
never caught earlier because this script had never actually been invoked at indices >1001 until
now). Deployed via `scp` (this file lives in the repo, not just on bluehive).

Resubmitted all 380 Mspec work units as one job, individually, generously sized:
`sbatch --partition=preempt --array=0-379%15 --time=03:00:00 --mem=120G --cpus-per-task=1
--export=ALL,IDX_OFFSET=760 submit_plain.sbatch` -> job 31344741. Verified the offset mechanism
directly once the first task ran: local index 0 + `IDX_OFFSET=760` -> global index 760 ->
`XVKIRI_XVMAGY__Mspec` -- correctly resolved to the *Mspec* technique, not single-taper, ruling out
an off-by-offset corruption risk before letting the rest run unattended. (That first task's log
read `skip (already done)` -- not a bug: `run_plain.py`'s own idempotency guard, most likely
tripped by a preemption/requeue cycle on `preempt`'s `PreemptMode=REQUEUE`, exactly the documented
behavior of that guard.)

Also hit, and want to flag rather than silently absorb: partway through this investigation,
`squeue` itself returned `slurm_receive_msg: Socket timed out` for several minutes (affecting both
manual checks and the persistent background monitor, which had been declaring `ALL_JOBS_DONE` off
of a misread empty/error `squeue` response -- a false positive, caught by cross-checking
`results_so_far` against the known 1520 total before trusting it). Replaced that monitor with a
version that treats a `squeue` error as "skip this tick," not "zero jobs," and additionally
requires `results_so_far >= 1520` before ever declaring completion.

Job 31344741 (380 Mspec units) is running; not yet confirmed complete as of this log entry.

### 2026-09-02 -- Stage 4: the multiprocessing driver retired for the rest of this batch --
FastMspec's timeouts kept spreading, and single-taper (assumed safe) turned out to be the worst
offender

While Mspec's plain-driver fixup (job 31344741) ran cleanly in the background, a routine sanity
sweep across the other three still-running mp-driver jobs (`sacct -j <job> --format=State -X | sort
| uniq -c`) turned up a much bigger problem than expected:

| Job (technique) | Completed | Timed out | Failure rate |
|---|---|---|---|
| 31344561 (FastMspec) | 45 | **25** (was 8 at last count) | ~36% of attempted |
| 31344563 (MspecBestK) | 66 | 2 | ~3% of attempted |
| 31344560 (single-taper) | 1 | **17** | **~89% of attempted** |

**FastMspec had kept failing even at the 2-hour extension.** The `scontrol update TimeLimit=02:00:00`
applied earlier (successfully, for then-pending tasks) was not enough -- 17 *more* tasks (`21, 23,
26, 28, 31, 33, 36, 38, 41-45, 51, 55, 58, 67`) timed out at `02:00:xx`, not just `01:00:xx`. Chunk=3
under real contention plus this technique's already-documented runtime variance (up to 3270s for a
single pair) can apparently exceed even a doubled budget when an unlucky combination of slow pairs
lands in the same chunk.

**Single-taper -- assumed the safe, fast case (2.47s in the Stage 3 pilot) -- was actually the
worst.** Its one successful task's own log (`mp_31344560_14.out`) shows why: real per-pair total
cost (cross-spectrum + the 33-template picking scan) ranged **186.6s to 984.8s**, not the 2.47s
pilot number -- that pilot only ever timed the raw cross-spectrum step in isolation, never the
picking/template-scan step that turns out to dominate total cost for every technique, single-taper
included. With `chunk=20`/`cpus-per-task=20` (real 20-way contention on one node) and a `00:30:00`
budget sized off the wrong (cross-spectrum-only) number, 17 of 19 tasks blew through it.

**Conclusion, not just a bigger patch:** the multiprocessing driver's core assumption --
same-node contention across `chunk` concurrent units, one shared time budget -- is unsafe for
*every* technique in this dataset, not only the expensive ones, because per-pair picking cost has
real, hard-to-predict variance no pilot sample fully captured. The plain driver (one work unit,
one dedicated task, one time budget) is structurally immune to this specific failure mode: a slow
pair costs only its own task's budget, never a neighbor's.

**Action taken**: retired the mp driver for all remaining/lost work across all four techniques,
funneling everything through `submit_plain.sbatch` from here on:
- Cancelled FastMspec's remaining pending mp tasks (`scancel 31344561_66 31344561_[73-126]`,
  before they could repeat the same failure) and let its 2 still-running tasks finish naturally.
- Cancelled MspecBestK's remaining pending mp tasks (`scancel 31344563_[73-126]`) pre-emptively --
  its failure rate was still low (~3%), but FastMspec looked fine early too (2 known losses) before
  ballooning to 25, so the same caution was applied rather than waiting for it to degrade the same
  way.
- Computed each technique's true missing-work-unit set directly (`build_work_units()` + checking
  `results/` for `work_unit_id + '.json'` -- **caught and fixed a bug in this exact check along the
  way**: an earlier draft of the script reconstructed the result filename from `pair.stn1`/
  `pair.stn2` directly, which lack the network-code prefix (`KIRI` vs. the real `XVKIRI`) that
  `work_unit_id` already carries correctly -- silently returned "380/380 missing" for FastMspec
  when 198 genuinely existed. Fixed by using `work_unit_id` directly; re-verified against a raw
  `ls results/*__FastMspec.json | wc -l` count (198) before trusting the corrected numbers.
- Resubmitted each technique's real missing set as one `submit_plain.sbatch` job, sized per
  technique from what's now known about its real cost:

| Technique | Missing | Job | `--time` | Notes |
|---|---|---|---|---|
| FastMspec | 182 | 31345341 | 02:00:00 | indices <1001, no offset needed |
| single-taper | 186 | 31345343 | 01:00:00 | indices <1001, no offset needed |
| MspecBestK | 168 | 31345366 | 01:00:00 | `IDX_OFFSET=1140` (indices exceed `MaxArraySize=1001`) |
| Mspec | (already running) | 31344741 | 03:00:00 | `IDX_OFFSET=760`, 0 failures so far |

All four now converging toward completion via the plain driver; none yet confirmed fully done as
of this log entry. **Follow-up for Stage 5**: the picking/template-scan step, not the cross-spectrum
computation, dominates real per-pair cost across every technique -- worth stating plainly in the
notebook's methods description rather than leaving the misleading impression (from the Stage 3
pilot table) that cross-spectrum cost alone characterizes technique expense.

### 2026-09-02 -- Stage 4: splitting remaining work across `preempt` and `standard`

Every job this whole batch has submitted so far ran on `preempt` (plus the one-off `debug` smoke
test) -- `standard` was never tried. Prompted by a direct question about this, checked its current
state: 91 nodes, mostly `mixed`/`allocated` but genuinely in use, 383 tasks queued from all users
combined and **zero from this account**. By contrast `preempt` had become substantially self-
congested -- at one point 688 tasks queued, 550 of them this account's own, competing against each
other under the same fair-share/fifo throttling. Per direct guidance, split the remaining pending
work in half for each of the four still-running plain-driver jobs: cancelled half of each job's
still-`PENDING` array indices on `preempt` and resubmitted that exact half as a fresh job on
`standard`, same per-technique `--time`/`--mem` budget as its `preempt` counterpart (no
`IDX_OFFSET` change needed -- same local-index/offset scheme, just a different `--partition`):

| Technique | preempt job (kept half) | standard job (moved half) | units moved |
|---|---|---|---|
| Mspec | 31344741 | 31345490 | 126 |
| FastMspec | 31345341 | 31345491 | 90 |
| single-taper | 31345343 | 31345492 | 93 |
| MspecBestK | 31345366 | 31345493 | 84 |

393 work units now running on two independent scheduling pools at once rather than one congested
one. `run_plain.py`'s own idempotency guard (`if out_path.exists(): skip`) means this split carries
zero risk of duplicate/wasted computation even though the same manifest indices are, in principle,
now known to two separate jobs.

### 2026-09-02 -- Stage 4: a first look at the partial results (742/1520, ~49%) -- three real findings

`aggregate.py` is explicitly designed to run on a partial batch (per its own docstring), so ran it
now rather than only tracking result *counts*. Real findings, not just a progress number:

**1. `single-taper`: 0/194 converged -- and it's a mechanistic zero, not just "worse."**
`freq_coverage_fraction` is **exactly 0.0 for all 193 non-null single-taper rows** (mean, std, min,
max all `0.0`) -- not low-and-variable, invariant. Cross-checked against the other diagnostics: 
single-taper has ~7-13x more `n_candidate_crossings` than any multitaper technique (1693.5 mean vs.
125-250) despite a *lower* `bad_quality_fraction` (0.74 vs. 0.87-0.93) -- consistent with an
unsmoothed, single-taper coherence estimate being far noisier (many more raw zero-crossings), which
then defeats the picker's sequential cycle-jump-tracking logic even though individual crossings
don't always trip the specific `bad_quality` gate. This is a real, physically-motivated result, not
an obvious bug -- and it's actually a strong, quantitative reinforcement of this whole project's
founding thesis (multitaper vs. single-taper), stronger than "somewhat better," closer to "the
picker structurally cannot extract a curve from single-taper's coherence at all" under the current
picker parameters. Flagging for Stage 5 rather than treating as settled: worth a plot of raw
crossing density (single-taper vs. multitaper) as direct visual evidence for this claim.

**2. `mean_amp_ratio` shows `inf` for 20-42% of multitaper work units, including converged ones --
a diagnostic-aggregation bug in this project's own Stage 2 instrumentation, not a correctness bug.**
Traced to `_vendored_seislib_an_processing.py:1008`:
```python
_diag_amp_ratios.append(maxamp / minamp if minamp > 0 else float('inf'))
```
Deliberate, not accidental -- when a pick's local envelope minimum is exactly zero the ratio really
is infinite, and that's physically plausible right near a Bessel/coherence null (exactly the region
this method targets). The bug is downstream: `_diag_mean_amp_ratio = np.mean(_diag_amp_ratios)`
means a *single* such pick poisons the whole work unit's mean to `inf`, discarding every other
finite, informative ratio from that unit. Confirmed this does **not** affect `converged`/
`best_delta_km_s` -- `_score()`'s `min(mean_amp_ratio / 5.0, 1.0)` saturates cleanly on `inf`, no
crash, no silent wrong answer in the actual picked results. It only makes `mean_amp_ratio` currently
unusable as a quantitative per-pair quality signal for a large fraction of units.
**Deliberately not patched mid-batch** -- doing so now would leave already-completed rows with the
old (broken) aggregation and new rows with a fixed one, an inconsistent column that would need
reconciling anyway. Recorded here as a concrete Stage 5 data-cleaning task instead: recompute a
robust central tendency from the raw per-pick ratios if they're ever needed (median instead of mean,
or finite-only mean + a separate reported null-pick-fraction), not a re-run of the batch.

**3. MATLAB cross-validation, now at real sample size (366/742 rows, vs. Stage 3's single pilot
pair) -- reassuring.** single-taper: mean relative L2 error **2.88e-10** (machine precision, 179
real pairs) -- confirms pipeline correctness at production scale, not just the one hand-picked
SKRH-BAND example. FastMspec: mean **6.38e-3** (0.6%), worst case 9.26% (`ANLA-ZOBE`) -- consistent
with, not new relative to, the already-documented and intentionally-not-reproduced MATLAB complex-
floor bug (`ccf_pipeline/NOTES.md`'s "Known upstream bug" section) characterized in Stage 3.

### 2026-09-02 -- Decision: two-round strategy, and what "saving everything" would need

Discussion prompted by a direct concern: this run's raw-data-to-coherence step is by far its most
expensive stage (memory and compute both), and right now only scalar diagnostics survive per work
unit -- if that step needs redoing to answer a question raised after the fact, the expensive part
is redone for nothing. Evaluated (not yet implemented) what a properly rich-data-capturing run
would need:

- **The real waste isn't the missing arrays, it's the template scan.** `work_unit.py`'s `process()`
  calls `extract_dispcurve()` up to ~33 times per work unit (once per corridor template,
  `n_templates_scanned` averages 32.8) and keeps only the best-scoring template's scalar summary --
  the other ~32 full picker runs vanish entirely, not even their scalars.
- Within even the kept run: the raw coherence curve (`faxis_pos`/`coh_pos`) is already in hand
  before the picker is even called (zero extra compute to save it); `zero_crossings` and
  `bad_quality`/`crossamps`/`peakamps` are computed unconditionally by `get_zero_crossings()`
  regardless of convergence, but only ever returned on the success path -- meaning failing units
  (the *majority*, 72-100% depending on technique) currently reveal nothing about why they failed,
  despite the informative arrays already existing in memory at failure time. All of this is a
  **zero-additional-compute, pure-I/O** fix.
- The one genuinely expensive item on the wishlist -- the KDE density grid + ellipse paths needed
  to regenerate seislib's own diagnostic plot -- only gets computed inside `extract_dispcurve`'s
  `if plotting:` block, which `work_unit.py` explicitly sets `False`. Capturing it is **real added
  compute**, not just an I/O change, and was recommended against as a blanket per-unit policy --
  better done on-demand for specific pairs of interest during Stage 5 (cheap at that scale).
- Full writeup, including size estimates (~75-230MB total for coherence+crossings+curves at the
  current 1520-unit scale, vs. 100MB-1.5GB for the KDE grid alone), given directly in conversation,
  not yet copied into a design doc.

**Decision**: let this run (Round 1) finish as designed -- not disrupted, its ~49%-and-climbing
progress preserved. Mine everything the *existing* scalar diagnostics can reveal in the meantime
(see below -- this already produced real findings). Once Round 1 completes, with full knowledge of
real per-technique timing distributions, the driver-design lessons (mp-driver retirement, the
`preempt`/`standard` split, `IDX_OFFSET`), and this diagnostics evaluation, design and run a
**Round 2**: a properly-sized, HPC-informed, rich-data-capturing full re-run, rather than patching
Round 1 mid-flight and paying for two half-informed runs instead of one well-designed one.
Rationale, in the user's own words: *"we are still in exploratory stage and compute runs are
expensive, so trying to understand why what we are doing works (or doesn't work) requires saving
almost everything. This is the key to debugging, and hypotheses testing."*

**Round 2 requirements (captured now, not yet designed in detail):**
- `WorkUnitResult` (or a richer sibling type) carries, per template attempt, not just the winner:
  `converged`, `bad_quality_fraction`, `freq_coverage_fraction`, `n_candidate_crossings`,
  `n_accepted_picks`, `best_delta_km_s` -- cheap, and answers the picking-stability-vs-reference-
  curve-choice question directly (see finding 2 below) instead of only via a converged-count proxy.
- Raw coherence curve, zero crossings (+ branch index), and `bad_quality`/`crossamps`/`peakamps`
  saved per (pair, technique) -- on both success and failure -- as a compact per-unit array file
  (`.npz`, not JSON) alongside the existing scalar JSON, matching the size estimates above.
- Driver: the decoupled persistent-worker/pull-queue architecture discussed for 1e5-1e6 scale
  (not required at 1520 units, but worth prototyping now while the problem is still small, given
  Round 2 is the natural place to validate it before it's needed for real).
- Consider widening the distance-quartile sweep beyond 4 bins, and adding a real
  convergence-rate-vs-distance figure to Round 2's own validation -- the resolution-bandwidth
  analysis below suggests it's worth extending, given how cleanly the theory already matches the
  Round 1 data.
- **Scope decision (2026-09-03), while Round 1 was still running): Round 2 tests only
  `single-taper` and `FastMspec`, dropping `Mspec`/`MspecBestK`.** Rationale, confirmed against
  Round 1's own data rather than assumed: FastMspec *is* Karnik/Romberg/Davenport's own algorithm
  (this project's actual reference implementation); Mspec/MspecBestK are Sayan's classical-
  multitaper comparison baselines, not a second "faithful" candidate. FastMspec is also the
  cheapest by construction (Mspec needed the K-chunking architectural fix mid-Round-1 and still
  drove most of this stage's memory/timeout bugs; MspecBestK's smaller-but-still-real K hit the
  same OOM class on the dataset's largest files, fixed twice this session). One refinement to the
  original reasoning ("all mspec techniques yield the same convergence rate"): they don't, quite --
  the distance-quartile table above shows Mspec notably underperforming FastMspec/MspecBestK at
  mid-range and MspecBestK slightly ahead at long-range -- doesn't change the choice (FastMspec
  still wins on both real criteria), just worth being precise about. Bonus reinforcement, from a
  different angle: `single-taper`/`FastMspec` are the *only* two techniques with real MATLAB
  cross-validation ground truth in this pipeline (`MATLAB_TECHNIQUE_DIRS` never covered Mspec/
  MspecBestK) -- this scope cut maximizes Round 2's validation power, not just its speed.

### 2026-09-02 -- Mining the existing scalars, per the two-round decision -- real findings without
any new data

**1. The resolution-bandwidth criterion (`2W ≲ Δf_zero = c/2r`) is directly confirmed in
real convergence data, at zero additional cost -- `dist_km` and `converged` were already in the
manifest.** Convergence rate drops monotonically with distance for every multitaper technique:

| Technique | Nearest quartile | 2nd | 3rd | Farthest quartile |
|---|---|---|---|---|
| FastMspec | 65% | 20% | 14% | 0% |
| Mspec | 60% | 6% | 0% | 0% |
| MspecBestK | 69% | 20% | 17% | 4% |

Exactly the predicted shape: `Δf_zero` shrinks as distance grows, so a fixed bandwidth `W`
that resolves a near pair's closely-spaced Bessel zeros becomes too coarse for a far pair. A real,
quantitative, dataset-derived validation of theory already established earlier in this project,
not a new claim -- worth a real figure in Stage 5's notebook (convergence-rate-vs-distance, one
line per technique).

**2. Template-scan stability, from the existing `n_templates_converged` count.** Among units where
the winning template converged, 40-55% of the ~33 other corridor templates also converged
(FastMspec 13.8/33, Mspec 18.1/33, MspecBestK 13.2/33) -- when a pair converges at all, it is not a
single-template fluke. Cannot yet say whether those other converged templates give *similar* or
*meaningfully different* `best_delta_km_s` values -- exactly the gap Round 2's per-template capture
(above) is meant to close.

**3. Converged units run ~15-20% longer** than non-converged ones, consistently across all three
multitaper techniques -- plausibly a richer-data proxy (more traces/longer duration -> both slower
and better coherence), not demonstrated as causal.

**4. MATLAB mismatch does not meaningfully predict non-convergence** (FastMspec: 0.10% median for
converged vs. 0.19% for non-converged, both small) -- and single-taper's non-converged units still
match MATLAB to 2.4e-12, reconfirming (again) that single-taper's 0% rate is a real methodological
result, not a computational error anywhere in the pipeline.

### 2026-09-02 -- Stage 4: one more Round-1 driver bug -- MspecBestK OOM on the dataset's largest
files; a local-machine session hiccup along the way

Picked back up after this machine's own Claude Code process restarted (the prior session's
background monitors died with it -- the bluehive jobs themselves kept running unattended the whole
time, unaffected). Briefly suspected a colliding parallel session when the tracker file showed
substantial content not in this restarted context's memory -- checked via `ListAgents`/
`SendMessage`: both other local sessions (`claude-sandbox-f7`, idle/unrelated; `claude-sandbox-e8`,
a different project, LitDiscovery, touching only its own `litdisc_judge_worker_*` jobs on
`urseismo`) confirmed no bluehive/tracker activity on this task. Root cause was simpler: this
restarted session's context predated the session's own most recent checkpoint, so the on-disk
tracker (and `git log`) were simply further along than this context's memory -- resolved by
re-reading current state and continuing from it, per direct guidance. One real mistake made along
the way and fixed immediately: an interrupted `Edit` call still partially wrote, corrupting the
tracker's own title line (`ost # Notebook 5 revamp...`) -- caught via `git diff` before committing
anything, reverted.

Resynced against real state: `sacct`'s per-job completion counts were reading anomalously low
across all four active jobs (not chased down -- possibly an accounting-DB windowing quirk on this
cluster) -- switched to counting `results/*__<technique>.json` directly instead, authoritative and
cheap. Current split: single-taper 209, FastMspec 216, Mspec 145, MspecBestK 218 -- 788/1520, ~52%.

A sweep of `.err` logs for crash signatures (`Bus error`, `Segmentation`, `Killed`, `OOM`) found one
live, new failure mode: **`31345366_83` and `31345366_175` (global indices 1223, 1315 --
`XVANLA_XVMMBE__MspecBestK` and `XVLONA_XVMAJA__MspecBestK`) both hit `Bus error (core dumped)` /
"Exceeded step memory limit"** against MspecBestK's `--mem=32G` budget. Their raw `.mat` files are
**1.5GB and 1.7GB** -- well above the dataset's ~690MB mean and bigger than anything Stage 3's own
timing/memory pilot (267MB SKRH-BAND) exercised, so `32G` was never actually validated against this
technique's real upper tail. Fixed the same way as every prior memory-sizing gap this stage:
resubmitted just these 2 indices with `--mem=120G` (already proven safe elsewhere) --
`sbatch --array=83,175 --time=01:00:00 --mem=120G --cpus-per-task=1 --export=ALL,IDX_OFFSET=1140
submit_plain.sbatch` -> job 31345552. A third crashed-looking log (`plain_31344555_760.err`,
`Killed`/OOM) turned out to be stale noise from an early single-pair Mspec smoke test (submitted
`01:09:52`, well before the real batch) whose pair (`XVKIRI_XVMAGY`) already has a genuine
completed result from a later run -- confirmed harmless, left alone.

Restarted persistent background monitoring (the prior monitor task IDs died with the session
restart) with one addition: each tick now also greps all `.err` logs for the same crash signatures
and reports a count, so a new OOM/segfault surfaces automatically rather than needing another
manual sweep to find.

### 2026-09-03 -- Stage 4: MspecBestK's `32G` budget was too small in general, not just for 2
outlier pairs

The crash-signature counter added above did its job: after another local session restart dropped
the background monitor again (bluehive jobs unaffected, as before), resyncing found the counter at
14, not the expected 3 -- 11 new crashes, all `Bus error`/OOM, all in job `31345366` (the original
MspecBestK plain-driver job, still at `--mem=32G`). So the fix for indices 1223/1315 two entries
ago was necessary but not sufficient -- `32G` was never actually a safe general budget for this
technique, just one that happened to survive the first ~39 attempts before hitting its next
large-file pair. `sacct -j 31345366` confirmed: 13 `FAILED` (the original 2 plus these 11) against
only 25 `COMPLETED` -- a ~34% failure rate on this job's attempted tasks, not a rare tail.

Fixed the same way as every prior mis-sizing this stage, and by now a routine move: cancelled every
remaining pending MspecBestK task at the old budget (`scancel 31345493_[296-379] 31345366_[250-
295]`, covering both the `preempt` and `standard` copies), recomputed the technique's true missing
set directly against `results/` (141 work units -- confirms the cancelled range was still mostly
unattempted, not wasted partial progress), and resubmitted all 141 in one job at the
already-proven-safe `--mem=120G`: `sbatch --array=<141 offset indices> --time=01:00:00 --mem=120G
--cpus-per-task=1 --export=ALL,IDX_OFFSET=1140 submit_plain.sbatch` -> job 31346905. No similar
crash signatures found for single-taper (still at `32G`) as of this check -- watching, not yet
assumed safe, given MspecBestK looked fine at the same budget for its own first ~39 attempts too.

### 2026-09-03 -- Why FastMspec doesn't need K-chunking, verified against the real code -- this
is the headline architectural argument for Stage 5, not just an empirical convenience

Prompted by a direct question: is "FastMspec never needed the K-chunking fix Mspec needed" a real,
citable architectural claim worth foregrounding in Stage 5, or just an empirical accident of this
project's parameter choices? Checked against the actual source rather than assumed.

**The mechanism** (`python/thomson_multitaper/fast_multitaper.py`, translated directly from
Santhosh Karnik's own `FastMultitaper.m`, 2019): classical `classical_spectrum_batch` (Mspec/
MspecBestK) does the obvious thing -- multiply data by all `K` tapers, FFT each, average. Memory
scales `O(N x n_traces x K)` by construction, exactly why K=80 pushed Mspec past 55GB before the
K-chunking fix. FastMspec instead splits the multitaper average into two pieces: `z0`, a *single*
FFT-based sinc-kernel convolution of the plain periodogram that analytically stands in for every
"well-concentrated" taper (DPSS eigenvalue ~1) at once, `O(N log N)` cost, never materialized
individually; and `z1`, an explicit correction over only the tapers in the DPSS eigenvalue
*transition region* (neither cleanly concentrated nor cleanly rejected) -- the code calls this
count `r` (`self.r = len(eig_weights)` in `fast_multitaper.py`), and it is a property of how wide
that transition band is, not of `K`.

**Confirmed `r` really does stay small/flat while `K` grows, rather than tracking it** (ran
`FastMultitaper` directly at N=10801 for several NW values):

| NW | K (full taper count) | r (transition-region correction) | r/K |
|---|---|---|---|
| 5 | 5 | 11 | 220% |
| 10.8 (this project's actual FastMspec bandwidth, `Wband=0.001`) | 15 | 13 | 87% |
| 20 | 33 | 14 | 42% |
| 50 | 92 | 16 | 17% |
| 100 (~ Mspec's own `K=80` setting) | 191 | 18 | **9%** |
| 200 | 390 | 20 | 5% |
| 400 | 790 | 20 | 3% |

Extended afterward (2026-09-03) across a wider NW=5-400 sweep: `r` stays essentially bounded
(**11-20**) across nearly two orders of magnitude of `NW`, while `K` grows linearly the whole way
(5 to 790). There's a crossover around `NW~10-20` -- coincidentally close to where this project's
current ad-hoc bandwidth already sits -- below which `r` can exceed `K` (the correction dominates,
no memory savings yet, though `r` itself is still trivially small in absolute terms), and above
which `r/K` collapses toward zero. This is the basis for the "coupled to the bandwidth search"
argument recorded in the deferred/requirement log above: wherever Stage 5's principled `W` search
lands, FastMspec's memory story is already settled.

**The honest nuance, stated explicitly rather than glossed over**: at this project's own narrow
bandwidth, `r` isn't dramatically smaller than `K` (13 vs 15) -- part of why Mspec looks so much
worse here isn't purely algorithmic, it's also that Mspec is configured with a much larger `K=80`
(`NW=100`) as a deliberately heavier classical baseline, not matched to FastMspec's resolution-
optimal bandwidth. Two effects, both real, worth separating in Stage 5:
1. **Architectural** (general, provably true at any bandwidth, the one to feature): `r` grows far
   slower than `K` -- the 9% ratio at NW=100 shows this isn't a narrow-bandwidth artifact, it holds
   *more* strongly as `K` grows, i.e. gets stronger exactly where memory would otherwise be worst.
   Traceable directly to Karnik et al.'s own method, not an empirical convenience of this project's
   choices.
2. **Configuration** (specific to this project's parameters, the honest caveat): Mspec's `K=80` is
   itself a larger taper count than FastMspec's own resolution-matched `Wband=0.001` implies, so
   the two techniques aren't only running different algorithms -- they're also targeting different
   taper counts. Keeps the architectural claim from overclaiming.

**Round 2 timing estimate, from Round 1's own real per-unit runtimes** (already representative --
solo, plain-driver, no shared-node contention):

| Technique | n sampled | mean | median | p90 | max | Full-380-pair total |
|---|---|---|---|---|---|---|
| single-taper | 239 | 474.8s | 457.1s | 758.2s | 1235.6s | ~50 CPU-hours |
| FastMspec | 258 | 1032.8s | 1033.1s | 1507.1s | 3269.9s | ~109 CPU-hours |

~159 CPU-hours total for Round 2's 760 work units (380 pairs x 2 techniques). Wall-clock depends
entirely on achievable concurrency, which Round 1 showed varies wildly on this shared cluster (2
to 30+ real concurrent tasks depending on fair-share/partition congestion) -- at a realistic
10-20 average, that's roughly 8-16 hours, with real risk of stretching past a day if the cluster
stays as congested as Round 1's worst stretches. Plan for "well under a day if it goes well, budget
up to ~2 days," not a single number.
