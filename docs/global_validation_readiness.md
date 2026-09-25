# Global FastMSPEC validation — readiness gates (implications of the pipeline smoke test)

*Written 2026-09-25. Program goal (Tolu): validate FastMSPEC at scale — ~2,000 stations worldwide (wavenet-epicAI `metadata3/fps_stations.csv`, fixed, never regenerate), global **Love-wave** noise correlations, **20 s dispersion maps compared with Ekström (GDM52)**, and a demonstration of FastMSPEC's global noise-correlation improvements.*

Source of the findings below: the XD.MTAN–XD.RUNG smoke test of the acquisition pipeline (`URseismology/wavenet-epicAI`, branch `add-september-ncf-pipeline`). The test lives in a **separate project**, `projects/wavenet_xd_pair_test/` (living handoff: `STATE.md`), and — once pushed — in `docs/ncf_pipeline_stages/xd_mtan_rung_smoke_test/` on that branch (report, patches, tests, evidence). This file records what it means **for FastMSPEC**.

## 1. Which existing FastMSPEC results are affected?
| result | data path | affected by the pipeline's day-placement defect? |
|---|---|---|
| Round 1 / Round 2 (380 `XV` pairs), Notebook 4 Round-1/2 sections | Sayan's matched data (`*_win_3_all_matched_data.mat`, MATLAB path) | **No** — different path |
| SKRH–BAND, BITY–MAGY benchmarks (Notebooks 3–4, `verification/nb3_technique_cost/`) | `ADAMA_gvib.h5` via `dispcurve_pick/gvib_loader.py` | **No.** Checked 2026-09-25 on AF.SKRH, XV.BAND, XV.BITY, XV.MAGY (60 day chunks each): each day is an independent array starting at an integer second (fractional start ≤ 0.025 s), so inter-station timing error ≈ 0.03 s versus ≈ 1 s in the orchestrator's continuous arrays. Note: some days carry 86,401 samples (the same "extra boundary sample" seen in the orchestrator data); `gvib_loader` anchors the last window to the shorter array's end, so both stations stay index-aligned. |
| Notebook 3 MTAN/RUNG (LH, SAC) | Sayan's SAC files | **No** |
| **Any global run fed by the orchestrator's `master.h5`** | continuous arrays, day-by-day appended | **Yes, unless the patches are adopted** |

The global set can only come through the orchestrator path (ADAMA_gvib covers ~1,372 stations, mostly Africa and permanent networks), so the pipeline fixes are on the **critical path** of the global program.

## 2. Gates before a global run (ordered by how much each can change a conclusion)
1. **Window length versus distance (independent of the pipeline defects).** The resolution ceiling `NW_high = N·c_min/(4R)` with the `K = 0` floor at NW ≈ 3 (Round 2, H3/H4). For illustrative c = 3.5 km/s:

   | R (km) | 1 000 | 2 000 | 5 000 | 10 000 | 15 000 | 20 000 |
   |---|---|---|---|---|---|---|
   | NW_high, N = 10 801 (3 h) | 9.5 | 4.7 | **1.9** | **0.9** | 0.6 | 0.5 |
   | NW_high, N = 86 401 (1 day) | 75.6 | 37.8 | 15.1 | 7.6 | 5.0 | 3.8 |

   Three-hour windows fall below the numerical floor beyond ~2 000 km — the same mechanism that gave quartile 4 (780–1 030 km) 0 % convergence in Round 2. A global run needs day-scale (or longer) windows, chosen per pair. This is what Notebook 05 (bandwidth selection, ADAMA benchmark) is meant to establish, and Round 2 §6's window-length recommendation (one splice test was negative — not yet a validated fix).
2. **The per-day edge taper becomes a first-order problem for long windows.** The orchestrator tapers 5 % at both ends of *every processed day* (72 min each). In a continuous master array that is a periodic 1/day amplitude notch covering ~10 % of the record; a day-scale window contains it in full. The smoke-test patches deliberately do **not** fix this (needs overlap-padded per-day processing, or longer continuous chunks). Gate: fix before day-scale windows are used.
3. **Adopt the pipeline fixes** (`xd_mtan_rung_smoke_test/PATCHES.md`): integer-second day alignment (removes a station-specific delay of k×0.05 s, up to ~1 s, that re-randomises after gaps), no dropped days (a day file with one extra sample made the whole next day be refused), response-derived units, per-day quality flags. At 20 s on long paths the timing error matters less than on this 110 km pair, but it is systematic, and **BH-derived stations carry it while native-1 Hz (LH) stations largely do not** — an unpatched network would be heterogeneous.
4. **Love-wave path is not yet validated.** The smoke test used BHZ (Rayleigh). Love needs the horizontals rotated to transverse with correct **sensor orientation**: `data/metadata/orientation.csv` already lists offsets for this very pair (XD-MTAN +2.47°, XD-RUNG −3.94°), whereas the orchestrator stores only StationXML's nominal azimuth. A 4° rotation error leaks sin 4° ≈ 7 % of the radial component into the transverse channel. Needs an orientation estimate for the 2 000 stations. **Cheap next test:** the packaged XD data already contain BHN/BHE; rotate → transverse coherence → benchmark against ADAMA `co_love` / `cf_love` for `XD.RUNG-XD.MTAN` (references already extracted in the smoke-test project).
5. **Quality control before stacking.** Corrupt raw days degrade coherence (8-day dry run: correlation with ADAMA's predicted coherence 0.49 → 0.56 after removing 3 of 8 days). A FastMSPEC-versus-conventional comparison must use identically screened data.
6. **Timing audit on a sample of the 2 000 stations** from raw headers before committing compute (`code/timing_replay.py` in the smoke-test folder needs only the day files' headers).
7. **Comparison target.** Our own `data/reference/hybrid_curve_README.md` records that GDM52 starts at **25 s**. Before promising a 20 s comparison with Ekström, confirm which model/period is available; otherwise compare at 25 s.
8. **Reference curve outside Africa.** The picker's per-pair hybrid reference curve takes ADAMA maps (6–40 s, **Africa only**) below 45 s and GDM52 from 45 s up (GDM52 itself starts at 25 s, but the hybrid leaves those periods to ADAMA). For a global pair outside Africa the 6–40 s part has no source, and 20 s lies below GDM52's shortest period, so the template corridor at 20 s would have to be extrapolated. This is the "global short-period reference-curve source" item deferred in the tracker (`docs/notebook5_revamp_progress.md`, "Deferred"), and it becomes a gate for the global run: extend `hybrid_reference_curve.py` to use GDM52 from 25 s where ADAMA is absent and extrapolate (or adopt a global short-period model), or run the picker with a wide, model-independent corridor and report picking success separately. The Rayleigh path is also still unvalidated against ADAMA.

## 3. What this does not change
FastMSPEC's per-pair findings (Notebooks 3–4, Round 1/2, Stage 4.5) stand: they were computed on data paths with no timing defect. The smoke test also **strengthens** one argument — FastMSPEC's cost stays flat as NW grows (Notebook 3 §1c) — which matters precisely because day-scale windows at NW 5–75 will be needed.

## 4. Suggested order
1. Land the pipeline fixes (push in progress; authors decide adoption). 2. Love-wave (transverse) smoke test on the same pair. 3. Timing audit on a station sample. 4. Notebook 05 → window/bandwidth design per distance. 5. Taper/overlap-padding fix. 6. Pilot (~50–100 stations), then the 2 000-station run.
