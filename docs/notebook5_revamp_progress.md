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
- [x] **Stage 1** -- Provenance/citations (pull 5 papers, `docs/references/README.md`, update
      `notebooks/README.md` + `docs/coherence_barcode_design.tex`)
- [x] **Stage 2** -- Vendor + instrument the seislib picker (`python/dispcurve_pick/`)
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

Next: Stage 3 (validate against Sayan's own SKRH-BAND result; 4-technique timing pilot).
