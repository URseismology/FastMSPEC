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
