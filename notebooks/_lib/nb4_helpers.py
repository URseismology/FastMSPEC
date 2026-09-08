"""Helper code for Notebook 5 (the coherence barcode). Full design rationale, citations, and
derivations: docs/coherence_barcode_design.tex / .pdf -- read that first if this module's
docstrings aren't enough context on their own.

Two halves:

1. Template barcodes (analytically exact): `load_reference_curve`, `build_template_family`, and
   `template_barcode` implement the Hawkins & Sambridge (2019) result that zero-crossings and
   extrema of J0(2*pi*f*r/c(f)) occur exactly at the tabulated zeros of J0 and J1 respectively --
   solved by root-finding against the *target* Bessel-zero value (exact), not by densely sampling
   the theoretical curve and numerically detecting sign changes (approximate).

2. Candidate barcodes (noise-contaminated) and matching: `candidate_barcode` extracts a
   reliability-filtered Z/M/N barcode from a real (noisy) coherency array, generalizing Notebook
   3 Section 1f's zero-crossing/local-amplitude-swing logic to maxima and minima. `match_score`
   and `scan_templates` implement the precision+recall+corridor-rejection+frequency-weighted
   scoring described in the design document.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.special import j0, jn_zeros


# ---------------------------------------------------------------------------
# Reference curve and template family
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Template barcodes: analytically exact (Hawkins & Sambridge 2019, Eqs. 1-3)
# ---------------------------------------------------------------------------

def _solve_bessel_events(c_interp: interp1d, dist_km: float, f_lo: float, f_hi: float,
                          zero_targets: np.ndarray) -> np.ndarray:
    """For a candidate phase-velocity curve c(f), solves x(f) = 2*pi*f*r/c(f) = target for f,
    for each target in zero_targets that actually falls within the curve's achieved x-range over
    [f_lo, f_hi]. x(f) is monotonic over any physically sensible dispersion curve in this band, so
    a single brentq bracket spanning the whole band suffices per target.
    """
    def x_of_f(f):
        return 2 * np.pi * f * dist_km / c_interp(f)

    x_lo, x_hi = x_of_f(f_lo), x_of_f(f_hi)
    lo, hi = min(x_lo, x_hi), max(x_lo, x_hi)
    events = []
    for target in zero_targets:
        if lo < target < hi:
            events.append(brentq(lambda f: x_of_f(f) - target, f_lo, f_hi, xtol=1e-12))
    return np.array(events)


def template_barcode(c_interp: interp1d, dist_km: float, f_lo: float, f_hi: float,
                      n_zeros: int = 500):
    """Analytically exact Z/M/N barcode for one template curve c(f) (Section 2.1 of the design
    doc). Returns a dict with 'Z' (zero-crossing frequencies), 'M' (maxima), 'N' (minima) -- M/N
    classified by evaluating J0 at each solved extremum against its immediate neighbors.
    """
    j0_zeros = jn_zeros(0, n_zeros)
    j1_zeros = jn_zeros(1, n_zeros)
    zc = _solve_bessel_events(c_interp, dist_km, f_lo, f_hi, j0_zeros)
    ext = _solve_bessel_events(c_interp, dist_km, f_lo, f_hi, j1_zeros)

    def x_of_f(f):
        return 2 * np.pi * f * dist_km / c_interp(f)

    eps = 1e-6 * (f_hi - f_lo)
    maxima, minima = [], []
    for f in ext:
        lo_val = j0(x_of_f(max(f - eps, f_lo)))
        hi_val = j0(x_of_f(min(f + eps, f_hi)))
        mid_val = j0(x_of_f(f))
        (maxima if (mid_val > lo_val and mid_val > hi_val) else minima).append(f)

    return {"Z": np.sort(zc), "M": np.sort(np.array(maxima)), "N": np.sort(np.array(minima))}


# ---------------------------------------------------------------------------
# Candidate barcodes: noise-contaminated, reliability-filtered
# ---------------------------------------------------------------------------

def _local_swing(y: np.ndarray, idx: np.ndarray, half_win: int = 3) -> np.ndarray:
    """Amplitude swing (max-min) in a small window straddling each event index -- generalizes
    directly across event types (Notebook 3 Section 1f used this for zero-crossings only)."""
    out = []
    for i in idx:
        lo, hi = max(0, i - half_win), min(len(y), i + half_win + 2)
        out.append(y[lo:hi].max() - y[lo:hi].min())
    return np.array(out)


def _zero_crossing_idx_and_freq(y: np.ndarray, freqs: np.ndarray):
    sgn = np.sign(y)
    idx = np.where(np.diff(sgn) != 0)[0]
    fc = freqs[idx] + (freqs[idx + 1] - freqs[idx]) * (0 - y[idx]) / (y[idx + 1] - y[idx])
    return idx, fc


def _extrema_idx_and_freq(y: np.ndarray, freqs: np.ndarray):
    """Local maxima/minima on a sampled (noisy) curve -- simple neighbor comparison, no
    closed form available for real data."""
    is_max = (y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])
    is_min = (y[1:-1] < y[:-2]) & (y[1:-1] < y[2:])
    max_idx = np.where(is_max)[0] + 1
    min_idx = np.where(is_min)[0] + 1
    return max_idx, freqs[max_idx], min_idx, freqs[min_idx]


def candidate_barcode(freqs: np.ndarray, coherency_real: np.ndarray, reliable_pct: float = 50.0):
    """Reliability-filtered Z/M/N barcode from a real (noisy) coherency array, using each event
    type's own within-method swing distribution as the reliability threshold (Section 2.2 of the
    design doc -- a shared absolute threshold across methods/event-types is meaningless, since
    e.g. raw single-taper's noise alone produces larger swings than FastMspec's entire signal).

    Returns a dict with 'Z'/'M'/'N' -> (all_freqs, reliable_mask) so callers can inspect both the
    full and reliability-filtered event sets.
    """
    zc_idx, zc_f = _zero_crossing_idx_and_freq(coherency_real, freqs)
    max_idx, max_f, min_idx, min_f = _extrema_idx_and_freq(coherency_real, freqs)

    zc_swing = _local_swing(coherency_real, zc_idx)
    max_swing = _local_swing(coherency_real, max_idx)
    min_swing = _local_swing(coherency_real, min_idx)

    def reliable_mask(swing):
        if len(swing) == 0:
            return np.array([], dtype=bool)
        thresh = np.percentile(swing, 100 - reliable_pct)
        return swing >= thresh

    return {
        "Z": (zc_f, reliable_mask(zc_swing)),
        "M": (max_f, reliable_mask(max_swing)),
        "N": (min_f, reliable_mask(min_swing)),
    }


# ---------------------------------------------------------------------------
# Matching and scoring
# ---------------------------------------------------------------------------

def _implied_velocity(f: float, event_freqs_sorted: np.ndarray, dist_km: float,
                       zero_index_table: np.ndarray) -> float | None:
    """Ekstrom et al. (2009)-style implied local velocity c = 2*pi*f*r/target, using the target
    Bessel-zero value nearest this event's own rank among event_freqs_sorted (a crude proxy for
    "which zero index is this" without needing to solve the full inverse problem)."""
    rank = np.searchsorted(event_freqs_sorted, f)
    if rank >= len(zero_index_table):
        return None
    target = zero_index_table[rank]
    if target <= 0:
        return None
    return 2 * np.pi * f * dist_km / target


def _one_to_one_match(cand_f: np.ndarray, template_f: np.ndarray, tol: np.ndarray):
    """Greedy one-to-one bipartite matching: each candidate event can match at most one template
    event and vice versa (nearest pairs assigned first). This is essential, not cosmetic -- a
    naive "does ANY candidate event fall within tolerance of this template event" check (and its
    mirror for precision) lets a dense candidate place multiple events near the same template
    event, inflating recall almost for free since precision doesn't penalize the redundancy. This
    was caught empirically during implementation: without one-to-one matching, raw single-taper's
    9x-denser event set outscored FastMspec outright, the opposite of every other result in this
    project. `tol` is per-template-event (already includes any weighting/scaling upstream).

    Returns (matched_template_idx, matched_candidate_idx) as boolean masks.
    """
    matched_t = np.zeros(len(template_f), dtype=bool)
    matched_c = np.zeros(len(cand_f), dtype=bool)
    if len(cand_f) == 0 or len(template_f) == 0:
        return matched_t, matched_c

    pairs = []
    for i, tf in enumerate(template_f):
        diffs = np.abs(cand_f - tf)
        within = np.where(diffs <= tol[i])[0]
        for j in within:
            pairs.append((diffs[j], i, j))
    pairs.sort(key=lambda p: p[0])

    used_t, used_c = set(), set()
    for _, i, j in pairs:
        if i not in used_t and j not in used_c:
            used_t.add(i)
            used_c.add(j)
    matched_t[list(used_t)] = True
    matched_c[list(used_c)] = True
    return matched_t, matched_c


def match_score(candidate: dict, template: dict, dist_km: float,
                 corridor_c_min: float, corridor_c_max: float,
                 tol_frac: float = 0.25, low_freq_weight: float = 2.0,
                 freq_split: float | None = None):
    """Precision+recall match score between a candidate barcode (from `candidate_barcode`,
    reliability-filtered) and one template barcode (from `template_barcode`), following the
    four-part design in Section 3 of the design document:

    1. Matching tolerance scaled to the template's own local event spacing (tol_frac of the gap
       to the nearest neighboring template event of the same type), not a fixed Hz value.
    2. One-to-one matching (see `_one_to_one_match`), then combined precision + recall (harmonic
       mean) over that matching -- not recall alone, and not a many-to-many "any nearby event
       counts" check either.
    3. Hard corridor rejection: a candidate event whose implied local velocity
       (Ekstrom et al. 2009-style c = omega*r/z) falls outside [corridor_c_min, corridor_c_max]
       is excluded from matching entirely, not merely counted as unmatched.
    4. Frequency-weighted trust: mismatches below `freq_split` (default: band midpoint) count
       `low_freq_weight` times as much as mismatches above it, following Hawkins & Sambridge's
       (2019) argument that reference-model deviation is "confined to the near surface or higher
       frequency."

    Returns a dict with the combined score, precision, recall.
    """
    total_weight = 0.0
    matched_weight = 0.0
    candidate_total_weight = 0.0
    candidate_matched_weight = 0.0

    for etype in ("Z", "M", "N"):
        template_f = template[etype]
        cand_f, cand_reliable = candidate[etype]
        cand_f = cand_f[cand_reliable]
        if len(template_f) == 0:
            continue

        f_lo_band = freq_split if freq_split is not None else 0.5 * (template_f.min() + template_f.max())

        sorted_t = np.sort(template_f)
        spacing = np.zeros(len(sorted_t))
        for i in range(len(sorted_t)):
            gaps = []
            if i > 0:
                gaps.append(sorted_t[i] - sorted_t[i - 1])
            if i < len(sorted_t) - 1:
                gaps.append(sorted_t[i + 1] - sorted_t[i])
            spacing[i] = np.mean(gaps) if gaps else np.inf
        tol = tol_frac * spacing

        zero_table = jn_zeros(0 if etype == "Z" else 1, max(len(sorted_t) + 5, 50))
        cand_kept = np.array([
            f for f in cand_f
            if (civ := _implied_velocity(f, sorted_t, dist_km, zero_table)) is not None
            and corridor_c_min <= civ <= corridor_c_max
        ])

        matched_t, matched_c = _one_to_one_match(cand_kept, sorted_t, tol)

        for f_t, is_matched in zip(sorted_t, matched_t):
            w = low_freq_weight if f_t < f_lo_band else 1.0
            total_weight += w
            matched_weight += w if is_matched else 0.0

        for f_c, is_matched in zip(cand_kept, matched_c):
            w = low_freq_weight if f_c < f_lo_band else 1.0
            candidate_total_weight += w
            candidate_matched_weight += w if is_matched else 0.0

    recall = matched_weight / total_weight if total_weight > 0 else 0.0
    precision = candidate_matched_weight / candidate_total_weight if candidate_total_weight > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {"score": f1, "precision": precision, "recall": recall}


def scan_templates(candidate: dict, template_library: dict, dist_km: float, f_lo: float, f_hi: float,
                    corridor_c_min: float, corridor_c_max: float, **score_kwargs):
    """Scores `candidate` against every template in `template_library` ({delta: c_interp}),
    returning (best_delta, best_score_dict, all_scores) where all_scores is {delta: score_dict}.
    """
    all_scores = {}
    for delta, c_interp in template_library.items():
        tpl_barcode = template_barcode(c_interp, dist_km, f_lo, f_hi)
        all_scores[delta] = match_score(candidate, tpl_barcode, dist_km,
                                         corridor_c_min, corridor_c_max, **score_kwargs)
    best_delta = max(all_scores, key=lambda d: all_scores[d]["score"])
    return best_delta, all_scores[best_delta], all_scores
