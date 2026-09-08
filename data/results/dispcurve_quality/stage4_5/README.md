# Stage 4.5 validation results

Two small (4-row) result tables from the modest-scale validation pass that closed out Stage 4.5
(`docs/notebook5_revamp_progress.md`'s 2026-09-07 log), both running `work_unit.process()`
(FastMspec, module-default `WBAND=0.001`, no `wband_override`) on the same 4 report example pairs,
one per distance quartile — isolating the hybrid-reference-curve fix's own effect, holding
bandwidth, corridor strategy, and the `horizontal_polarization` fix fixed across both runs.

- **`stage4_5_validation_old_curve.csv`** — bluehive job `31351550`: the single shared `SDISPL.ASC`
  reference curve (pre-Stage-4.5 default), combined with the corridor-widen strategy and
  `horizontal_polarization` fix. Reconstructed from that job's own stdout log
  (`stage4_5_validation_31351550.out`, printed via `df.to_string()`) since its CSV output path was
  later overwritten by the job below (both wrote to the same filename on bluehive) -- not a
  re-derivation, a direct parse of the job's own recorded output.
- **`stage4_5_validation_new_curve.csv`** — bluehive job `31351552`: the per-pair hybrid
  ADAMA+GDM52 curve (`hybrid_reference_curve.py`), same corridor/polarization fixes. Pulled
  directly as the CSV `work_unit.py`'s own driver script (`validate_stage4_5.py`) wrote.
- **`corridor_strategy_eval.csv`** — bluehive job `31351542`: the earlier 3-way corridor-strategy
  comparison (BASELINE / widen-at-short-period / down-weight-in-scoring) that chose the
  widen-at-short-period strategy used in both files above. Same 4 pairs, uses the hybrid curve
  throughout (this comparison predates and is independent of the old-vs-new-curve isolation above).

Headline result (see `docs/notebook5_revamp_progress.md`'s 2026-09-07 log for the full writeup):
the number of corridor templates that converge at all rises with the hybrid curve in every one of
the 4 pairs (0->8, 2->15, 15->25, 2->9 out of 33 scanned), with the winning template's own quality
(coverage, bad-quality fraction) also improving on 3 of 4 pairs.

Small enough to commit directly (~1KB each), unlike the raw matched-data inputs these came from.
