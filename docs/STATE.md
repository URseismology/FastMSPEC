# STATE — FastMSPEC status board (living; update at every milestone)

*Last updated: 2026-09-25. Branch `notebook5-phase-velocity-revamp`, pushed through `beb28f0`.
Entry point for agents: [`../CLAUDE.md`](../CLAUDE.md). Catalog of every document: [`INDEX.md`](INDEX.md).*

The project has several work streams that advanced independently. This board says, per stream, **where it stands, what is open, and where the detail lives**.
Status words: **done** (verified, committed), **awaiting review** (finished, waiting for Tolu's feedback), **open** (next work), **blocked** (needs something first), **deferred** (recorded, not scheduled).

## Program goal
Validate FastMSPEC at scale: ~2,000 stations worldwide → global **Love-wave** noise correlations → global **20 s dispersion maps compared with Ekström (GDM52)** → demonstrate FastMSPEC's improvement.
Readiness gates and the implications of the pipeline smoke test: [`global_validation_readiness.md`](global_validation_readiness.md).

## Streams

| # | stream | status | what is open | detail |
|---|---|---|---|---|
| A | **Multitaper library** (`python/thomson_multitaper`) | done | — | [`../python/NOTES.md`](../python/NOTES.md), [`../verification/octave_verify_multitaper/README.md`](../verification/octave_verify_multitaper/README.md) |
| B | **CCF pipeline translation** (`python/ccf_pipeline`) | done (six phases; Octave-verified) | one unresolved discrepancy in the Butterworth filter helper; memory fix `classical_spectrum_batch` done 2026-09-01 | [`../python/ccf_pipeline/NOTES.md`](../python/ccf_pipeline/NOTES.md), [`plan_ccf_mtc_translation.md`](plan_ccf_mtc_translation.md) |
| C | **Dispersion-curve revamp, Stages 0–5** (vendored picker, Round 1, Round 2, Stage 4.5 hybrid reference curve + corridor + Love polarization fix) | done | Stage 6 only: final pass on `notebooks/README.md` status column/timings | [`notebook5_revamp_progress.md`](notebook5_revamp_progress.md) ("Status summary" at top; dated log below) |
| D | **Round 1 / Round 2 results** (380 pairs × 4 techniques; 300-point bandwidth sweep) | done — **historical**, computed with the pre-Stage-4.5 picker | a re-run with the current picker was deliberately not scheduled | [`round2_hypothesis_evaluation.pdf`](round2_hypothesis_evaluation.pdf) (H1, H3, H4 confirmed; H2 not supported; §6.5: the FastMSPEC argument is architectural, not "always dramatically cheaper") |
| E | **Notebooks 01, 02** | done | — | [`../notebooks/README.md`](../notebooks/README.md) |
| F | **Notebook 03** `03_fastmspec_application` | **awaiting review** — Sections 1–2b and the Summary reworked and executed for real 2026-09-09 (commits `105d090`, `7dbc8c3`; 1,870 s, 0 errors, 6 figures): technique cost, memory, convergence and coherence quality on AF.SKRH–XV.BAND and XV.BITY–XV.MAGY, NW sweep anchored to the resolution ceiling, MTAN/RUNG demoted to 2b | Tolu's feedback on the executed notebook; Sections 2c (second pair, Bessel-fit metric), 2d (envelope conditioning) and 3 (synthetic NLNM) are carried over from the earlier version and were not reworked | source `../notebooks/_lib/build_nb3.py`; precompute drivers [`../verification/nb3_technique_cost/README.md`](../verification/nb3_technique_cost/README.md); artifacts `../data/results/dispcurve_quality/nb3_precompute/` |
| G | **Notebook 04** `04_dispersion_curve_picking` | done, executed for real (`cff0115`, 607 s, 13 figures) after the figure-feedback revision | — | tracker "RESOLVED 2026-09-08" block |
| H | **Notebook 05** `05_bandwidth_selection` (bandwidth-selection theory `NW_low` vs the ADAMA benchmark) | **open — not started** | follow "Start here" in the ADAMA tracker: (1) finish overlap check for the 141 island-touching pairs, (2) parameterize the windowing convention in `gvib_loader.build_pair_matched_data`, (3) run all 141 pairs, (4) FastMspec vs ADAMA `co`/`cf`, (5) window length N versus distance, (6) build `notebooks/_lib/build_nb5.py` | [`findLowBand_ADAMAbenchmark_progress.md`](findLowBand_ADAMAbenchmark_progress.md), [`stage5_bandwidth_theory.pdf`](stage5_bandwidth_theory.pdf) |
| I | **Notebook 06** coda correlation | scaffold only; deliberately out of scope | — | [`../notebooks/README.md`](../notebooks/README.md) |
| J | **Global validation program** (2,000 stations) | **blocked** on the gates below | Love-wave (transverse) smoke test; timing audit on a station sample; window/bandwidth design per distance (Notebook 05); taper/overlap padding; reference curve outside Africa; pilot of 50–100 stations | [`global_validation_readiness.md`](global_validation_readiness.md) |
| K | **Acquisition-pipeline smoke test** (WaveNet-EpicAI, XD.MTAN–XD.RUNG) — *adjacent project, separate folder* | fixes verified on a 4-day window and 4 hard cases; **full-window runs in progress on bluehive**; report/PDF/push to `URseismology/wavenet-epicAI` branch `add-september-ncf-pipeline` pending | pull results → figures → `REPORT.md`/`REPORT.pdf` → assemble sub-folder → push (authorized by Tolu) | local `~/claude-sandbox/projects/wavenet_xd_pair_test/STATE.md` (not in this repository) |

## Implications of stream K for streams F–J (why the separate folder still matters here)
* **Existing FastMSPEC results are not affected.** Round 1/2 and the SKRH–BAND / BITY–MAGY benchmarks run on Sayan's matched data or on `ADAMA_gvib.h5` (each day an independent array starting at an integer second, fractional start ≤ 0.025 s); the smoke-test defect is in the orchestrator's continuous, day-appended arrays. Details: [`global_validation_readiness.md`](global_validation_readiness.md) §1.
* **The global run can only use the orchestrator path**, so the pipeline fixes (integer-second alignment, no dropped days, response-derived units, quality flags) are on its critical path.
* **Window length is the biggest scientific gate**, independent of the pipeline: 3-hour windows cannot resolve beyond ~2,000 km; day-scale windows make the per-day 5 % taper a first-order problem. This is the design question Notebook 05 is meant to answer.
* **Love waves are unvalidated end to end** (orientation, transverse rotation, picker on new data); Rayleigh is unvalidated against ADAMA in the picker.

## Decisions that should not be re-litigated
| decision | source |
|---|---|
| `ADAMA_gvib.h5` is the data source for ADAMA benchmarks; SAC files are already response-corrected | [`findLowBand_ADAMAbenchmark_progress.md`](findLowBand_ADAMAbenchmark_progress.md) |
| Notebooks use real data, best/hardest examples, show-before-tell, seislib-style KDE picking figures | Tolu, 2026-09-08/09 (Notebook 3 rework) |
| Use a precomputed solution plus commented-out reproduce commands when a classical Mspec run is slow | Tolu, 2026-09-09 |
| Notebook 3 must test large NW (e.g. 50) **against the ceiling** `NW_high` from `c_min`, `r`, theory | Tolu, 2026-09-09 |
| Do not claim "no major bugs" for the acquisition pipeline; report the timing defect | 2026-09-25 |
| Notebook numbering (01–06) and the `notebook5-v1-event-scanning` tag | tracker, 2026-09-08 |

## Deferred (recorded so nothing is silently lost)
* M/N (maxima/minima) quality gate for a three-way barcode (only zero-crossings are used).
* A global short-period reference-curve source (see gate 8 in the readiness document).
* Rayleigh-wave validation of the hybrid reference curve.
* `gvib_loader.chunk_days` memory fallback (raises `NotImplementedError`).
* Overlap-padded per-day taper in the acquisition pipeline; logger/inspector completeness checks (stream K recommendations).
* Notebook 03 Sections 2c, 2d and 3 (see stream F): re-evaluate against the reworked story once Tolu has reviewed 1–2b.

## Housekeeping (as of this update)
* `CLAUDE.md`, `docs/STATE.md`, `docs/INDEX.md`, `docs/check_links.py` and `docs/global_validation_readiness.md` were pushed to the remote in `beb28f0` (2026-09-25, Tolu's OK).
* **Decision (Tolu, 2026-09-25): no folder reorganization**; the entry-point documents are the fix. Options that were considered are at the end of [`INDEX.md`](INDEX.md).
* **Notebook 03 review is pending; Tolu is busy elsewhere and will come back to it.** Do not rework Sections 2c, 2d or 3 until that feedback arrives.
* The tracker's "Status summary" is dated 2026-09-08 and predates the Notebook 3 rework and stream K; this board supersedes it as the cross-stream entry point, and the tracker remains authoritative for streams C and D.

## How to update this file
One row per stream; change the status cell and add date + commit id; if a stream has its own tracker, add the dated entry there and keep this row to one sentence.
Run `python3 docs/check_links.py` afterwards.
