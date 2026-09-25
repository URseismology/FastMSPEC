# Documentation index

Agent entry point: [`../CLAUDE.md`](../CLAUDE.md). Cross-stream status: [`STATE.md`](STATE.md). This file catalogs every document and says how much to trust it.

**Status key.** *current* = maintained, safe to act on · *living* = updated every milestone · *historical* = a true record of what was done at the time, not the current pipeline · *reference* = background (papers, design rationale) · *stale-header* = content is valid but the top-of-file status text predates later work.

## Top level
| document | status | what it is |
|---|---|---|
| [`../CLAUDE.md`](../CLAUDE.md) | current | agent entry point: layout, rules, environments, git conventions |
| [`../README.md`](../README.md) | current (see note) | project overview and where-to-start reading order; its status paragraph points to the tracker and to `STATE.md` |
| [`STATE.md`](STATE.md) | living | cross-stream status board (what is done / open / blocked) |
| [`global_validation_readiness.md`](global_validation_readiness.md) | current (2026-09-25) | gates for the 2,000-station global run and what the pipeline smoke test implies for FastMSPEC |

## Trackers and plans (`docs/`)
| document | status | what it is |
|---|---|---|
| [`notebook5_revamp_progress.md`](notebook5_revamp_progress.md) | living log, *stale-header* (dated 2026-09-08) | dispersion-curve revamp, Stages 0–6, Round 1/2 and Stage 4.5; **2,285 lines** — read the "Status summary" at the top, search the dated log for reasons |
| [`findLowBand_ADAMAbenchmark_progress.md`](findLowBand_ADAMAbenchmark_progress.md) | current | plan, findings and "Start here" for Notebook 05 (bandwidth selection vs the ADAMA benchmark); data-source decision (`ADAMA_gvib.h5`) |
| [`plan_ccf_mtc_translation.md`](plan_ccf_mtc_translation.md) | historical (all six phases done 2026-08-18) | why the translation targets `ccf_compute_crosscorr_mtc_Z/T.m`, phased plan |
| [`manuscript_intro_draft.md`](manuscript_intro_draft.md) | reference (running draft) | opening paragraph of the eventual manuscript |

## Theory and reports (`docs/*.tex`, `.pdf`)
| document | status | what it is |
|---|---|---|
| [`round2_hypothesis_evaluation.pdf`](round2_hypothesis_evaluation.pdf) (`.tex`) | historical results (pre-Stage-4.5 picker); conclusions cited widely | Round 2 bandwidth-sweep evaluation of hypotheses H1–H4; §5.1 resolution ceiling `NW_high`; §6.5 architectural argument for FastMSPEC |
| [`stage5_bandwidth_theory.pdf`](stage5_bandwidth_theory.pdf) (`.tex`) | reference; **theory not yet checked against real data** (that is Notebook 05's job) | `NW_low` bias–variance bandwidth-selection theory |
| [`coherence_barcode_design.pdf`](coherence_barcode_design.pdf) (`.tex`) | reference; the original barcode-matching metric was replaced by the picker-based metric (its §8 documents the revision) | design of the zero-crossing/peak "barcode" and the instrumented reference-curve approach |
| [`references/`](references/README.md) | reference | the papers this work builds on (Karnik et al., Haley & Anitescu, Hawkins & Sambridge, Ekström, Magrini/seislib, Park et al., …) |
| [`figures/`](figures/) | historical | figures for the Round 2 report (`round2_report/`) and investigations (`investigation/`) |

## Code documentation
| document | covers |
|---|---|
| [`../python/NOTES.md`](../python/NOTES.md) | `thomson_multitaper`: translation notes, Octave verification, MEX audit |
| [`../python/ccf_pipeline/NOTES.md`](../python/ccf_pipeline/NOTES.md) | `ccf_pipeline`: status per phase, design decisions, upstream bug, memory fix |
| [`../python/dispcurve_pick/NOTES.md`](../python/dispcurve_pick/NOTES.md) | vendored/instrumented seislib picker: provenance, instrumentation, environment fix |
| [`../python/dispcurve_pick_batch/NOTES.md`](../python/dispcurve_pick_batch/NOTES.md) | bluehive batch pipeline (380 pairs × 4 techniques, Round 2 sweep) |
| [`../notebooks/README.md`](../notebooks/README.md) | the six notebooks, status, how to run them |
| [`../data/README.md`](../data/README.md), [`../data/reference/README.md`](../data/reference/README.md), [`../data/reference/hybrid_curve_README.md`](../data/reference/hybrid_curve_README.md), [`../data/results/dispcurve_quality/README.md`](../data/results/dispcurve_quality/README.md) | data provenance, reference curves, batch result layout |

## Verification (`../verification/`)
| folder | verifies |
|---|---|
| [`octave_verify_multitaper/`](../verification/octave_verify_multitaper/README.md) | Python multitaper vs the unmodified MATLAB source, run in Octave |
| [`octave_verify_ccf_pipeline/`](../verification/octave_verify_ccf_pipeline/README.md) | Python CCF pipeline vs the unmodified MATLAB pipeline, incl. real-data test |
| [`skrh_band_real_data/`](../verification/skrh_band_real_data/README.md) | first large real-data end-to-end check (SKRH–BAND, N = 10,801, 1,605 traces) |
| [`gvib_skrh_band_test/`](../verification/gvib_skrh_band_test/README.md) | `ADAMA_gvib.h5` as a data source (rotation bug and zero-day dilution bug found and fixed) |
| [`jackknife_variance_check/`](../verification/jackknife_variance_check/README.md) | closed-form expected jackknife variance (Haley & Anitescu) |
| [`window_length_splice_test/`](../verification/window_length_splice_test/README.md) | splice test for longer windows (negative result) + Round 2 report figures |
| [`nb3_technique_cost/`](../verification/nb3_technique_cost/README.md) | precompute drivers behind Notebook 3 Sections 1–2 (cost, memory, convergence, coherence) |

## Legacy (read-only history)
[`../legacy/matlab_source/README.md`](../legacy/matlab_source/README.md) (which `.m` file backs which Python module), [`../legacy/para_ccf_original/NOTE.md`](../legacy/para_ccf_original/NOTE.md) (the repository's pre-translation contents).

## Outside this repository (adjacent work)
| item | where | notes |
|---|---|---|
| Acquisition-pipeline smoke test (XD.MTAN–XD.RUNG) | `~/claude-sandbox/projects/wavenet_xd_pair_test/` (own git repository), `STATE.md` inside | will also be published to `URseismology/wavenet-epicAI`, branch `add-september-ncf-pipeline`, under `docs/ncf_pipeline_stages/xd_mtan_rung_smoke_test/` |
| Bluehive batch root | `/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch` | Round 1/2 results, Notebook 3 precompute |
| Bluehive smoke-test root | `/scratch/tolugboj_lab/wavenet_ncf_xd_pair_test/` | never touch `wavenet_ncf_production` / `wavenet_ncf_canary` |

## Proposed reorganization (not executed — needs Tolu's decision)
Measured cost of moving anything: the tracker is referenced from **27** other files (notebooks, `NOTES.md` files, README), `round2_hypothesis_evaluation` from 13, `stage5_bandwidth_theory` from 10, `coherence_barcode_design` from 9, `plan_ccf_mtc_translation` from 8, `findLowBand_ADAMAbenchmark_progress` from 5. Executed notebooks embed some of these paths in saved output text.

1. **Recommended — no file moves.** The navigation problem was missing entry points, not misplaced files. `CLAUDE.md`, `STATE.md`, `INDEX.md` and `check_links.py` now provide them at no cost to existing links.
2. **Optional, if you want a cleaner `docs/`:** group into `docs/trackers/` (the two trackers, plan), `docs/reports/` (the three `.tex`/`.pdf` pairs), leaving `STATE.md`, `INDEX.md`, `global_validation_readiness.md`, `references/`, `figures/` in place. Requires rewriting about 70 references and re-running `check_links.py`; also breaks any external link to the current paths (the GitHub repository may be linked from elsewhere).
3. **Optional, split the tracker:** freeze `notebook5_revamp_progress.md`'s dated log (lines after "Log") into an archive file, keeping the Status summary, Checklist and Deferred log as the living part. Same reference cost as moving the file unless the living part keeps the current name (then only the log moves and references to specific dates would need updating).
4. **Housekeeping to do regardless:** refresh the root `README.md` status paragraph and the tracker's "Status summary" header to point at `STATE.md` (done additively in this change), and update `notebooks/README.md`'s Status column for Notebook 03 after Tolu's review.
