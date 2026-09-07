"""Reference-curve loading and the AkiNet-corridor template family.

Moved here (Stage 4 of the Notebook 5 revamp) from `notebooks/_lib/nb5_helpers.py`, which
originally implemented these against the old zero-crossing/max-min event-scanning barcode design
-- unchanged logic, just relocated so both the bluehive batch pipeline (`dispcurve_pick_batch`)
and the (Stage 5, rewritten) notebook can import from one place, without the batch pipeline
depending on `notebooks/_lib` (not importable standalone on bluehive). Full design rationale:
docs/coherence_barcode_design.tex.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def load_reference_curve(path: Path, freqmin: float, freqmax: float, lmode: int = 0):
    """Parses SDISPL.ASC (Herrmann CPS-format multi-mode dispersion ASCII), filters to a single
    mode (fundamental, lmode=0, by default -- matching Sayan Swar's own use of this file), and
    restricts to [freqmin, freqmax]. Returns an interp1d over that restricted band, plus the
    band's actual (data-native) min/max frequency -- these can differ slightly from freqmin/
    freqmax themselves, since the file's own frequency grid rarely lands exactly on round numbers,
    and any solver using this curve must stay within the interpolator's actual domain.
    """
    df = pd.read_csv(path, sep=r"\s+", header=0)
    df.columns = [c.strip().replace("(", "_").replace(")", "").replace("/", "_") for c in df.columns]
    mode = df[df["LMODE"] == lmode].sort_values("FREQUENCY_Hz")
    band = mode[(mode["FREQUENCY_Hz"] >= freqmin) & (mode["FREQUENCY_Hz"] <= freqmax)]
    if len(band) < 2:
        raise ValueError(f"Fewer than 2 points for LMODE={lmode} in [{freqmin}, {freqmax}] Hz")
    c_interp = interp1d(band["FREQUENCY_Hz"], band["C_KM_S"], kind="linear")
    return c_interp, float(band["FREQUENCY_Hz"].min()), float(band["FREQUENCY_Hz"].max())


def build_template_family(c_ref: interp1d, f_lo: float, f_hi: float,
                           corridor_km_s: float = 0.8, step_km_s: float = 0.05):
    """c_template(f) = c_ref(f) + delta, for delta spanning +/-corridor_km_s in step_km_s steps
    (default: AkiNet's own tuned +/-0.8 km/s corridor, Xue & Olugboji 2025 Section 3.3.1 --
    adopted directly, not re-derived). Returns a dict {delta: interp1d} -- each value callable
    the same way as c_ref itself, over the same [f_lo, f_hi] domain.
    """
    f_grid = np.linspace(f_lo, f_hi, 500)
    c_ref_vals = c_ref(f_grid)
    deltas = np.arange(-corridor_km_s, corridor_km_s + step_km_s / 2, step_km_s)
    return {float(d): interp1d(f_grid, c_ref_vals + d, kind="linear") for d in deltas}


def build_template_family_widened(c_ref: interp1d, f_lo: float, f_hi: float,
                                    corridor_km_s: float = 0.8, step_km_s: float = 0.05,
                                    caution_period_s: float = 12.0, short_multiplier: float = 3.0,
                                    floor_km_s: float = 0.3):
    """"Widen the corridor at low-confidence (short) periods" -- Stage 4.5's corridor-strategy
    evaluation (docs/notebook5_revamp_progress.md, 2026-09 log; scratchpad
    eval_corridor_strategies.py, bluehive job 31351542) tested this against the unmodified
    build_template_family (BASELINE) and a scoring-side down-weight alternative on the 4 report
    example pairs. Result: this strategy was never worse than baseline (Q1/Q2: same diagnostics,
    Q4: still non-converged either way -- that quartile's failure is a separate, deeper issue, not
    a corridor-width problem, consistent with the branch-continuity/reseed negative results in
    docs/round2_hypothesis_evaluation.tex) and was a clear win on Q3 (coverage 0.904->0.923,
    bad_quality 0.125->0.062 vs. baseline). The down-weight alternative was never better than this
    one where they differed. Adopted as the corridor strategy used by work_unit.py.

    Same delta range/step/template count as build_template_family (identical search cost) -- only
    each template's low-confidence (period = 1/f < caution_period_s) portion gets
    `delta * short_multiplier` instead of `delta`, a proportionally wider excursion at short period
    for the same search cost. Physical floor (floor_km_s) clips the result: for a low reference
    curve combined with a large negative delta*short_multiplier, c_ref + delta*mult can go
    negative -- an unphysical phase velocity that was confirmed (via bluehive job 31349488's 7-hour
    zero-output timeout, diagnosed with a per-template timing script) to make extract_dispcurve
    hang indefinitely rather than fail fast. Clipping to a small positive floor keeps the
    "widen at low-confidence periods" concept intact while ruling out that pathological case.
    """
    f_grid = np.linspace(f_lo, f_hi, 500)
    c_ref_vals = c_ref(f_grid)
    short_mask = (1.0 / f_grid) < caution_period_s  # low-confidence = short period
    multiplier = np.where(short_mask, short_multiplier, 1.0)
    deltas = np.arange(-corridor_km_s, corridor_km_s + step_km_s / 2, step_km_s)
    return {float(d): interp1d(f_grid, np.clip(c_ref_vals + d * multiplier, floor_km_s, None),
                                kind="linear") for d in deltas}
