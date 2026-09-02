# Vendoring notes: `seislib.an.extract_dispcurve` -> `dispcurve_pick`

[← Back to repo README](../../README.md) | See also: [docs/coherence_barcode_design.tex](../../docs/coherence_barcode_design.tex) (Section 8, "Revision: From Barcode Matching to Instrumented Reference-Guided Picking")

Same role as [`ccf_pipeline/NOTES.md`](../ccf_pipeline/NOTES.md) for that package: tracks exactly
how the vendored source differs from upstream, and why.

## Provenance

`_vendored_seislib_an_processing.py` is Sayan Swar's own copy of seislib's internal picking
module, pulled from the lab NAS (`repovibranium`) at
`/volume1/web/FastMSPEC_data/madagascar_codes/python/_an_processing_seislib_functions.py`, via
`ssh repovibranium "cat <path>"` (same restricted-SFTP workaround used for `SDISPL.ASC` and the
`docs/references/` papers). 1539 lines.

Diffed directly against the installed `seislib==1.2.1` package's own
`seislib/an/_an_processing.py` at pull time: identical except one line (`np.in1d` vs. `np.isin` --
see "Environment-compatibility fix" below) and a trailing newline. Sayan's copy is, for all
practical purposes, exactly upstream seislib's own source, not a modified fork of it.

## What's instrumented, and why

`seislib.an.extract_dispcurve` computes rich internal quality signals while picking a dispersion
curve -- a per-crossing `bad_quality` flag (peak-ratio, envelope-relative-amplitude,
spacing-vs-reference criteria), the local kernel-density amplitude ratio behind every accepted
pick, and a frequency-coverage fraction its own acceptance test already checks against a 1/5
threshold -- but discards all of it. The public return is binary: `(crossings, dispersion_curve)`
on success, or a bare `DispersionCurveException` on failure, with no way for a caller to
distinguish *why* picking failed or *how good* a successful pick really was.

Three sites are instrumented, each wrapped in a `# --- FastMSPEC instrumentation: begin/end ---`
comment block, none of which alter picking behavior:

1. **`bad_quality` capture** (right after the 3-criterion gate + propagation-smoothing loop
   finishes): stashes `bad_quality`'s mean (the bad-quality fraction) and `len(w_axis)` (candidate
   crossing count) before `bad_quality` is consulted throughout the rest of picking below.
2. **Per-pick amplitude ratio capture** (inside the nested `pick_velocity` closure, immediately
   after `picks.append([frequency, vpick, maxamp])`): appends `maxamp/minamp` -- exactly the ratio
   the `pick_threshold` gate a few lines above just tested this pick against -- to an
   outer-scope list, via the same closure-mutation pattern the original code already uses for
   `picks` itself (no `nonlocal` needed, since only `.append()` is called, never a rebind).
3. **The final return/raise site**: adds an opt-in `return_diagnostics=False` parameter to
   `extract_dispcurve`'s signature (last positional/keyword slot, so no existing call site's
   argument order breaks). When `False` (the default), behavior is byte-identical to unmodified
   upstream -- same 2-tuple return, same plain `DispersionCurveException` on failure (see
   `tests/test_matches_upstream.py`). When `True`, returns a 3-tuple
   `(crossings, dispersion_curve, PickDiagnostics(...))` on success, or raises
   `DispersionCurveExceptionWithDiagnostics` (a `DispersionCurveException` subclass carrying a
   `.diagnostics` attribute) on failure -- deliberately not silently discarding diagnostics on the
   failure path, since a large sweep is expected to have many non-converging pairs, and knowing
   *why* (crossings were bad-quality vs. simply too few of them) is itself useful signal.

`freq_coverage_fraction` reuses the exact `len(smooth_picks) / len(logspaced_frequencies)`
quantity the original acceptance test already computes and compares against `1/5` -- not
recomputed independently. `PickDiagnostics` and `DispersionCurveExceptionWithDiagnostics` live in
`diagnostics.py`, imported into the vendored file rather than defined inline, to keep the vendored
file's own diff against upstream as small and legible as possible.

**One known gap**: `manual_picking=True` (an interactive `ginput`-based mode) returns early via an
un-instrumented code path -- calling with `manual_picking=True, return_diagnostics=True` together
silently returns the plain 2-tuple, not a 3-tuple. Not fixed, since manual picking is irrelevant to
this project's batch use case and is never exercised with `return_diagnostics=True` anywhere in
this repo -- flagged rather than silently accepted.

## Why the `seislib` package dependency was removed (Stage 4)

`_vendored_seislib_an_processing.py` originally imported 3 utility functions and 3 exception
classes from the installed `seislib` package (`from seislib.utils import ...` /
`from seislib.exceptions import ...`). Discovered during Stage 4 (deploying to bluehive): `pip
install seislib` there fails to build a completely unrelated Cython extension
(`seislib.tomography._ray_theory._math`, pre-generated `.c` code with a genuine variable-
redeclaration conflict against this toolchain's gcc) -- `pip` compiles every extension in a
package as part of building its wheel, so this blocks installing `seislib` at all, even though
dispersion-curve picking never touches `seislib.tomography`.

Fix: vendor the 3 functions (`adapt_timespan`, `adapt_sampling_rate`, `running_mean`, plus
`resample`, a helper `adapt_sampling_rate` itself calls) and 3 exception classes
(`DispersionCurveException`, `TimeSpanException`, `NonFiniteDataException`) directly --
`_vendored_seislib_utils.py`, `_vendored_seislib_exceptions.py`. Confirmed functionally identical
to the installed package via a direct source diff at vendoring time (one function,
`DispersionCurveException`, is literally byte-identical; the rest differ only in cosmetic operator
spacing and some trimmed docstring "Notes" sections -- no logic changed). This is a genuine
dependency-surface reduction, not just a bluehive workaround: the picking path no longer needs the
full `seislib` package (with its unrelated tomography/plotting submodules) installed at all,
anywhere. `seislib` itself is still a real, used dependency elsewhere in this repo (Notebook 3
Section 4's direct `seislib.an` import, and `tests/test_matches_upstream.py`'s own byte-fidelity
check against real upstream) -- only this module's production picking path was decoupled from it.

`tests/test_matches_upstream.py` was updated accordingly: it now asserts our own vendored
`DispersionCurveException` (not literally `seislib.exceptions.DispersionCurveException`, a
distinct-but-functionally-identical class since we no longer import it) is raised with an
identical message to what upstream would raise, checked by string comparison rather than
`isinstance`.

## Environment-compatibility fix (not a behavior change)

Sayan's copy uses `np.in1d` (line ~796, inside the manual-picking branch's crossing-lookup logic),
which numpy 2.x has fully removed (`AttributeError: module 'numpy' has no attribute 'in1d'`,
confirmed on this environment's numpy 2.4.6). The currently-installed `seislib==1.2.1` package
itself already uses `np.isin` at the equivalent line -- the same operation, just the
non-deprecated numpy 2.x spelling -- so this is not a translation deviation from real upstream
behavior, only from the specific (older-numpy-era) snapshot Sayan's own copy happened to be pulled
from. Fixed by swapping to `np.isin`, matching what upstream itself already does.

## Relationship to Notebook 3 Section 4's existing `seislib` usage

Notebook 3 Section 4 (`notebooks/_lib/build_nb3.py`) imports and calls
`seislib.an.extract_dispcurve` directly from the installed pip package -- unmodified upstream, no
instrumentation, no access to any of the diagnostics this package adds. That call site is
unaffected by anything here and is not in scope for this revamp (per the design doc's own stated
boundary). Notebook 5 uses this package (`dispcurve_pick`) instead, specifically for the
diagnostics `seislib.an.extract_dispcurve` alone cannot provide.

## Verification

`tests/test_matches_upstream.py` (3/3 passing): a synthetic Bessel-model coherence spectrum
(`Re{coherence(f)} = J0(2*pi*f*r/c(f))`, no noise -- a clean case is exactly right for a
translation-fidelity check, as distinct from a real-data validation) run through both the
unmodified installed `seislib` package and this vendored+instrumented copy, asserting:
- `return_diagnostics=False` (default): byte-identical `(crossings, dispersion_curve)` output.
- `return_diagnostics=True`: crossings/curve still byte-identical; the new `PickDiagnostics` is
  well-formed (fractions in `[0, 1]`, positive counts/ratios) on a converging case.
- A deliberately-too-narrow frequency band: both the unmodified package and this package
  (`return_diagnostics=False`) raise the plain `DispersionCurveException`; this package with
  `return_diagnostics=True` raises `DispersionCurveExceptionWithDiagnostics` carrying
  `diagnostics.converged is False`.

A real-data validation -- reproducing Sayan Swar's own SKRH-BAND FastMspec-vs-first-order result
from `phasevel_compute_slide13and14.ipynb` -- is Stage 3 of the revamp plan, not this stage; it
also doubles as the empirical timing pilot for sizing the bluehive batch job.
