"""Thin helpers for Notebook 4, written fresh for the phase-velocity/instrumented-picker
framework -- not retrofitted from the old event-scanning barcode's helpers (which this file
replaces; the old Z/M/N template-matching scorer -- candidate_barcode, match_score,
_one_to_one_match, scan_templates, template_barcode, _solve_bessel_events, and the raw
zero-crossing/extrema extraction helpers -- is removed outright per the approved Stage 5 plan,
recoverable via the `notebook5-v1-event-scanning` git tag, not worth keeping dead alongside its
own documented insufficiency, docs/coherence_barcode_design.tex Section 8).

`load_reference_curve`/`build_template_family` also moved out of this file during Stage 4 --
they're canonical now in `python/dispcurve_pick/template_family.py`, imported from there directly
by anything that needs them (this notebook included), not duplicated here.

A custom 2-way good/bad-quality-crossing barcode (`plot_quality_barcode`) was built, used, then
removed the same day (2026-09-08), per direct user feedback on the built notebook: `extract_dispcurve`'s
own `plotting=True` native visualization (reference curve, tracked branch, kernel-density field,
low-quality crossings all in one purpose-built figure) is a real improvement over that invention,
not a style preference -- see the notebook's own Figures 4/5 for the replacement. Recoverable via
git history if ever needed again, not worth keeping dead alongside its own superseded rationale.
"""
from __future__ import annotations

import signal

import numpy as np

from dispcurve_pick import extract_dispcurve, DispersionCurveExceptionWithDiagnostics

# A single extract_dispcurve call given a pathological template (near-zero/negative velocity
# range) can hang far past any reasonable per-call cost -- confirmed directly during Stage 4.5
# (bluehive job 31349488 ran 7 hours with zero output before this exact cause was diagnosed; see
# docs/notebook5_revamp_progress.md's 2026-09-06 log and build_template_family_widened's own
# floor_km_s clip, which prevents the specific cause found there but isn't a substitute for this
# defense-in-depth timeout at the call site itself). Same pattern as
# eval_corridor_strategies.py/work_unit.py's own PER_CALL_TIMEOUT_S.
PER_CALL_TIMEOUT_S = 90


class _PickTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _PickTimeout()


def scan_templates_with_picker(faxis_pos, coh_pos, dist_km, templates, f_lo, f_hi,
                                 cmin, cmax, horizontal_polarization=True,
                                 filt_width=10, filt_height=1.0, x_step=0.05, pick_threshold=0,
                                 verbose=True):
    """Runs the instrumented picker once per template in `templates` ({delta: interp1d}, from
    build_template_family/build_template_family_widened), returning a list of
    (delta, curve_or_None, picks_or_None, diagnostics_or_None) -- diagnostics is None only if
    this specific call hit PER_CALL_TIMEOUT_S, otherwise always present (even on a non-converged
    attempt) since every call passes return_diagnostics=True. This is the same
    scan-every-template-then-pick-the-best-scoring-one pattern work_unit.py's process() uses in
    the real bluehive batch pipeline; here it exists so the notebook can run a small, live
    demonstration of that same logic on one real pair, without needing the batch pipeline itself.
    """
    out = []
    for delta, c_template in templates.items():
        freqs_grid = np.linspace(f_lo, f_hi, 200)
        ref_curve_arr = np.column_stack([freqs_grid, c_template(freqs_grid)])
        signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(PER_CALL_TIMEOUT_S)
        try:
            curve, picks, diag = extract_dispcurve(
                faxis_pos, coh_pos, dist_km, ref_curve_arr,
                freqmin=f_lo, freqmax=f_hi, cmin=cmin, cmax=cmax,
                filt_width=filt_width, filt_height=filt_height, x_step=x_step,
                pick_threshold=pick_threshold, horizontal_polarization=horizontal_polarization,
                manual_picking=False, plotting=False, return_diagnostics=True,
            )
        except DispersionCurveExceptionWithDiagnostics as e:
            curve, picks, diag = None, None, e.diagnostics
        except _PickTimeout:
            curve, picks, diag = None, None, None
        finally:
            signal.alarm(0)
        if verbose:
            conv = diag.converged if diag is not None else 'TIMEOUT'
            print(f"  delta={delta:+.2f}: converged={conv}")
        out.append((delta, curve, picks, diag))
    return out


def score(diag) -> float:
    """Same combined score work_unit.py's own `_score` uses -- convergence gates everything,
    then rewards frequency coverage, low bad-quality fraction, and a strong mean amplitude ratio.
    Duplicated here (not imported) since work_unit.py's `_score` is a private, undecorated
    module-level function, not part of dispcurve_pick_batch's public surface -- and the formula
    is simple enough that duplication is more honest than reaching into another package's
    internals for something this small.
    """
    if diag is None or not diag.converged:
        return 0.0
    bad_q_term = 1.0 - min(diag.bad_quality_fraction, 1.0)
    amp_term = min(diag.mean_amp_ratio / 5.0, 1.0)
    return diag.freq_coverage_fraction + 0.5 * bad_q_term + 0.5 * amp_term


def best_of(scanned):
    """Picks the highest-scoring (delta, curve, picks, diag) tuple from scan_templates_with_picker's
    output; returns (None, None, None, None) if nothing converged."""
    best = max(scanned, key=lambda t: score(t[3]))
    return best if score(best[3]) > 0 else (None, None, None, None)
