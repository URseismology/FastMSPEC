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

- [ ] **Stage 0** -- This tracking file
- [ ] **Stage 1** -- Provenance/citations (pull 5 papers, `docs/references/README.md`, update
      `notebooks/README.md` + `docs/coherence_barcode_design.tex`)
- [ ] **Stage 2** -- Vendor + instrument the seislib picker (`python/dispcurve_pick/`)
- [ ] **Stage 3** -- Validate against Sayan's SKRH-BAND result; 4-technique timing pilot
- [ ] **Stage 4** -- bluehive batch pipeline, full 380 pairs x 4 techniques
- [ ] **Stage 5** -- Notebook 5 complete overhaul (built fresh, old version tagged not deleted)
- [ ] **Stage 6** -- Packaging + docs cleanup, Notebook 3 Section 4 ref_curve fix

## Deferred / requirement log

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

Next: Stage 1 (provenance/citations).
