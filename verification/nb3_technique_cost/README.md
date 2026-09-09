# Notebook 3 Section 1-2 precompute: technique cost & coherence

Generates the small artifacts Notebook 3's Sections 1, 1b, 1c and 2a load
(`data/results/dispcurve_quality/nb3_precompute/`). The cross-spectra need the
full ~250 MB matched-data arrays and the slow classical-Mspec runs, so they are
computed once on bluehive; the notebook does all picking and plotting live on the
small 1-D coherence arrays.

Runs on bluehive (`/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch`), in the
`fastmspec_batch` conda env, against `ADAMA_gvib.h5`.

| script | sbatch | output | notebook use |
|---|---|---|---|
| `driver_sec1.py` | `submit_sec1.sbatch` | `nb3_sec1_skrh.npz`, `nb3_sec1_metrics.json` | Section 1 (4 techniques, production configs) + 1c (FastMspec vs MspecBestK cost sweep, NW 10.8-40) |
| `driver_sec1b.py` | `submit_sec1b.sbatch` | `nb3_sec1b_ksweep.json` | Section 1b (classical Mspec K-sweep at fixed NW: roughness + runtime + peak memory) |
| `driver_sec2.py` | `submit_sec2.sbatch` | `nb3_sec2_bity_magy.npz` | Section 2a (XV.BITY-XV.MAGY, single-taper vs FastMspec coherence) |
| `worker.py` | — | (called per run by the drivers) | one technique / one bandwidth in a fresh process, so `ru_maxrss` peak memory is a clean isolated number |

**Deploy + run:**

```bash
scp verification/nb3_technique_cost/*.py verification/nb3_technique_cost/*.sbatch \
    bluehive:/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/nb3_precompute/
ssh bluehive 'cd /scratch/tolugboj_lab/FastMSPEC_dispcurve_batch && \
    sbatch nb3_precompute/submit_sec1.sbatch && \
    sbatch nb3_precompute/submit_sec1b.sbatch && \
    sbatch nb3_precompute/submit_sec2.sbatch'
# then pull report_figures/nb3_sec1*.{npz,json} and nb3_sec2_bity_magy.npz
# into data/results/dispcurve_quality/nb3_precompute/
```

## Config notes (so the numbers match the batch pipeline)

- `Mspec` in Section 1 uses `NW=100 / K=80` -- `work_unit.py`'s `NW_MSPEC` / `K_MSPEC`,
  the deliberately-heavy classical baseline behind every manifest and Round 2 number.
  `MspecBestK` and `FastMspec` use `Wband=0.001` (the production bandwidth).
- Section 1c sweeps only the two bandwidth-tracking techniques (`FastMspec`,
  `MspecBestK`) over NW in {10.8, 18, 24, 28, 40}. 10.8-28 stays at or below
  SKRH-BAND's resolution ceiling `NW_high = N*c_min/(4R)` (~25-30 for a physical
  Love-wave `c_min` ~3 km/s); NW=40 is deliberately past it, labelled as a cost
  probe, not a usable estimate.
- Matched data is rebuilt from `ADAMA_gvib.h5` via `dispcurve_pick.gvib_loader`
  (the same validated window+rotate path as `verification/gvib_skrh_band_test/`),
  giving `coh_num` 1620 for SKRH-BAND vs. Stage 3's precomputed-`.mat` 1605 --
  the small difference is the zero-filled-day exclusion, expected.

## Result (2026-09-08 run)

SKRH-BAND (289 km, built from gvib, 108 days, coh_num 1620), production configs:

| technique | K | runtime | peak mem | coherence range | pickable |
|---|---|---|---|---|---|
| single-taper | 1 | 3.9 s | 1.7 GB | [-0.122, 0.107] | no |
| MspecBestK | 15 | 112 s | 10.9 GB | [-0.092, 0.085] | yes |
| FastMspec | 13 | 137 s | 17.2 GB | [-0.092, 0.083] | yes |
| Mspec (NW=100) | 80 | 564 s | 11.7 GB | [-0.021, 0.016] | no (washed flat) |

Section 1c NW sweep (runtime / K):

| NW | FastMspec | MspecBestK |
|---|---|---|
| 10.8 | 137 s / K=13 | 112 s / K=15 |
| 18 | 144 s / K=14 | 210 s / K=29 |
| 24 | 134 s / K=14 | 277 s / K=41 |
| 28 | 144 s / K=14 | 350 s / K=49 |
| 40 (past ceiling) | 165 s / K=16 | 509 s / K=72 |

FastMspec runtime flat (~135-165 s, r stays ~13-16); MspecBestK climbs 4.5x with K.
Memory: both bounded -- MspecBestK ~12 GB (K-axis chunking cap), FastMspec ~18 GB
(higher, from the O(N log N) sinc term, but flat by construction).
