"""Builds 04_coda_correlation_future_work.ipynb from scratch via nbformat.
This is a roadmap/scaffold notebook, not a full implementation -- see the
plan's Notebook 4 section for why."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text))


nb['cells'] = cells

md(r"""# Notebook 4 — Future Work: Coda-Correlation

Sayan's report is explicit that this was out of scope for the course
project: *"Due to time constraints the project has not been extended to
computing the Coda Correlograms... This will be addressed in future works."*
No MATLAB script or Python code for it exists anywhere in this codebase (one
partial, non-functional draft was found — see Section 3). This notebook is a
**roadmap, not an implementation**: it documents the methodology precisely
enough that picking it up later doesn't require re-deriving it from scratch.
""")

md(r"""## 1. Theoretical grounding: the coda-correlation wavefield

Source: Tkalčić, Phạm & Wang (2020), "The Earth's coda correlation
wavefield: Rise of the new paradigm and recent advances," *Earth-Science
Reviews* 208, 103285.

Ordinary ambient-noise cross-correlation (Notebooks 1-3's subject) relies on
a **diffuse ambient wavefield** — many uncorrelated, roughly
equally-distributed sources — for the cross-correlation between two
receivers to approximate the inter-receiver Green's function:

$$
[G(x_B, x_A; t) + G(x_B, x_A; -t)] * S_N(t) = f_A(t) \star f_B(t)
$$

(Tkalčić et al.'s Eq. 1) — this is exactly the physical basis Notebook 2's
Section 1 built the Bessel-coherence relationship on, and it is fundamentally
a **surface-wave** phenomenon at the scales ambient-noise studies typically
operate at.

**Earthquake coda** — the part of a seismogram recorded many hours after the
first body-wave arrivals — offers a different, complementary wavefield.
Tkalčić et al. describe the coda-correlation wavefield's defining feature as
containing **many cross-terms of reverberating body waves**, a principle
"fundamentally different from the reconstruction of surface waves in the
ambient-noise correlograms." Where ambient-noise CCFs are strongest for
imaging shallow structure via surface waves, coda cross-correlations open a
route to probing **deeper** structure (the paper's abstract specifically
motivates this by ambient noise's inability to reach the depths that
receiver-based/body-wave studies can) using the same cross-correlation
machinery this repo already implements — provided the input windows are
chosen from **post-earthquake coda**, not arbitrary continuous noise.

This is why coda-correlation is a natural *extension* of the FastMspec work
in Notebooks 1-3, not an unrelated project: the stability gains (Notebooks
1, 3) and computational tractability (Notebooks 1-2) that FastMspec brings
to ambient-noise cross-spectra apply equally to computing cross-spectra of
coda windows — the spectral-estimation math doesn't care what physical
regime the input windows come from.
""")

md(r"""## 2. Sayan's stated plan (transcribed precisely from the report)

From the "How do you implement in Code" and "Results and Discussions"
sections of `SSWAR_ESC425_Project_Report.pdf`:

1. **Event-conditioned windowing.** Filter noise windows to only those
   starting **at least 2 hours after** moderate-to-large earthquake events
   (**Mw > 5.5**) occurring during the network's operational period. This
   is a fundamentally different window-selection criterion than
   `prepare_data.py`'s existing day/window pairing, which has no concept of
   event timing at all — every window is currently treated identically
   regardless of what happened nearby in time.
2. **Two parallel stacks for direct comparison**, in the report's own words:
   *"I will probably keep two sets, one is with seismic events and one
   without. Will compare the NCF wavefield for both of them. Seismicity
   Stack and Quiet Stack comparison."*
3. **Full stack, not sub-windowed.** Explicitly *not* separate 3-6h/6-9h
   post-event sub-stacks — *"Do full stack. Not just 3-6 or 6-9 hours stack
   separately. Choose window after earthquake and do full stack of all
   those windows."*
4. Two secondary, lower-priority items also noted in the report: NCF
   computation with 1-degree spatial binning ("if possible"), and a
   slowness-histogram 2D spatial-direction plot for introducing the dataset.
""")

md(r"""## 3. What's genuinely missing to execute this

A full-tree search for any coda-specific processing code found exactly one
relevant file:
[`ccf_time_spec_normalization_codacorr.m`](../legacy/matlab_source) (found
in `sayan-swar-translation/.../codes/test/matlab/lib/`, not currently copied
into this repo's `legacy/` since it isn't part of the verified translation
target) — but on inspection, it is an **incomplete, non-functional draft**:
line 6 references an undefined variable (`dt` instead of `parameters.dt`),
and line 23 (`amp=abs(datafft);phi`) has a syntax error — a stray token left
mid-edit. It implements generic time/spectral normalization, not
event-conditioned window selection, and would not run as-is even for that
narrower purpose. This confirms rather than contradicts the plan's original
assessment: **no working coda-correlation code exists in this codebase**,
draft or otherwise.

Two concrete, currently-missing pieces, in the order they'd need to be built:

1. **An earthquake catalog cross-referenced against station-pair operational
   windows.** Needed: event origin times and magnitudes (Mw) for the
   Madagascar network's operational period, from a standard catalog (e.g.
   USGS ComCat / ISC), filtered to Mw > 5.5, then intersected with each
   station pair's actual day-coverage (which `prepare_data.find_day_pairs`
   already determines) to know which windows are even eligible.
2. **Event-proximity window-filtering logic**, distinct from
   `prepare_data.build_windows`'s existing sliding-window cut logic. This
   needs to: (a) take the catalog from step 1, (b) for each qualifying
   event, locate the first window starting ≥2 hours after origin time, (c)
   select *all* subsequent windows through the end of that day's available
   data (per the report's "full stack, not sub-windowed" instruction) as
   the "Seismicity Stack" input, and (d) mark everything else as "Quiet
   Stack" input for the parallel comparison.

Both are new logic, not present in `ccf_pipeline` in any form — building
them is future work, not something this documentation pass should improvise
without an actual event catalog and Sayan's confirmation of catalog source
and Mw threshold interpretation (e.g. does "at least 2 hours after" mean
strictly after the event's P-wave arrival at the *station*, or after origin
time globally? The report doesn't specify, and it matters for correctness).
""")

md(r"""## 4. Concrete next step

For whoever picks this up next (Sayan, the user, or a future session):

1. **Get an earthquake catalog.** The USGS ComCat API
   (`https://earthquake.usgs.gov/fdsnws/event/1/`) can be queried directly
   for Mw > 5.5 events within the Madagascar network's operational date
   range — this is a single scripted API call, not a data-engineering
   project, and would be the natural first step.
2. **Resolve the "2 hours after" ambiguity** with Sayan directly (origin
   time vs. estimated station arrival time) before writing the filtering
   logic, since it changes which windows land in which stack.
3. **Write `find_coda_windows(catalog, station_pair_windows, min_delay_hours=2)`**
   as a new function in `ccf_pipeline` (not this notebook), following the
   same translation-quality bar as the rest of this repo — since here there
   is no MATLAB source to translate or verify against, this would need its
   own test suite (synthetic catalog + synthetic windows, checked by hand)
   rather than an Octave comparison.
4. Once windows are correctly split into "Seismicity Stack" and "Quiet
   Stack" sets, **everything downstream is already built**: feed both sets
   through the exact same `compute_crosscorr_mtc_fastmspec` pipeline
   Notebook 3 already exercises on real data, and compare the two resulting
   NCF wavefields exactly as Notebook 3 Section 2 compares single-taper vs.
   FastMspec — same comparison pattern, different input-selection criterion.

This notebook stops here deliberately — implementing step 1 onward is a new
work item, not a continuation of the documentation pass Notebooks 1-3
completed.
""")

with open('04_coda_correlation_future_work.ipynb', 'w') as f:
    nbf.write(nb, f)
print("wrote 04_coda_correlation_future_work.ipynb with", len(cells), "cells")
