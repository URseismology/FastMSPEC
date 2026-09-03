"""Round 2 preparation: select the stratified-by-distance NW-sweep subset, and compute each
selected pair's own NW_high(R) reference from its real Round 1 FastMspec result.

Not part of the driver itself -- a one-off analysis script, run once against Round 1's completed
results to produce the subset manifest + fraction-grid work list Round 2's driver will consume.
See docs/notebook5_revamp_progress.md's 2026-09-03 entries for the full derivation and rationale:
- NW_high(R) ~ N*c_min/(4*R), c_min taken over the pair's own analysis band from its Round 1
  reference-curve-plus-best_delta_km_s reconstruction (the correct, dispersion-aware choice --
  NOT a single assumed-constant c, and NOT the band's mean/max c, which would be under-conservative).
- Fraction grid [1.5, 1.0, 0.7, 0.45, 0.25] x NW_high(R), applied per pair, not a shared global NW.
- Subset: 15 pairs per distance quartile (60 total), quartiles computed over the FULL 380-pair
  manifest (not just currently-completed pairs), matching Round 1's own convergence-vs-distance
  quartile method for direct comparability.

Known limitation, stated plainly rather than glossed over: Round 1 saves only a scalar
`best_delta_km_s` per pair, not the picker's own traced curve -- so the "c(f)" used here is the
*template* curve (`c_ref(f) + best_delta_km_s`) that scored best, not the raw picked curve's own
local wiggles. A reasonable first-order proxy, not a substitute for Round 2's own richer per-unit
curve capture. Likewise, "most stable band" here means "the pair's full analysis band" -- Round 1 has no
per-frequency reliability signal to identify a true stable sub-band; that granularity is a Round 2
capability (once return_diagnostics=True's per-pick data is saved), not available yet. Where a
pair has no converged Round 1 pick to source `c(f)` from (common, deliberately not excluded --
see `c_min_source` below), the bare reference curve (delta=0) is used as a physically-reasonable
fallback rather than dropping the pair from the candidate pool.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from dispcurve_pick import load_reference_curve

FRACTION_GRID = [1.5, 1.0, 0.7, 0.45, 0.25]
N_PER_QUARTILE = 15
N_SAMPLES = 10801  # this project's fixed window length (window3hr convention)
PICK_FREQMIN, PICK_FREQMAX = 0.01, 0.5


def compute_c_min(best_delta_km_s: float, ref_curve_path: Path) -> float:
    """c_min over [PICK_FREQMIN, PICK_FREQMAX] for the winning template curve
    c_ref(f) + best_delta_km_s -- the conservative (worst-point-in-band) reference
    NW_high(R) needs, not the band's mean or max.
    """
    c_ref, f_lo, f_hi = load_reference_curve(ref_curve_path, PICK_FREQMIN, PICK_FREQMAX)
    f_grid = np.linspace(f_lo, f_hi, 500)
    c_vals = c_ref(f_grid) + best_delta_km_s
    return float(c_vals.min())


def nw_high(c_min_km_s: float, dist_km: float) -> float:
    """NW_high(R) ~ N * c_min / (4*R). c in km/s, R in km -> c/R in Hz, consistent with this
    project's WBAND convention (dt=1.0s sampling, frequencies in Hz-equivalent).
    """
    return N_SAMPLES * c_min_km_s / (4 * dist_km)


def select_subset(manifest_csv: str, results_dir: str, ref_curve_path: str) -> pd.DataFrame:
    """Returns one row per selected pair: dist_km, quartile, best_delta_km_s, c_min, nw_high,
    and the 5 target NW values (fraction_grid x nw_high) Round 2 should run FastMspec at.
    """
    manifest = pd.read_csv(manifest_csv).rename(columns={"stndist": "dist_km"})
    manifest["quartile"] = pd.qcut(manifest["dist_km"], 4, labels=[1, 2, 3, 4])

    results_dir = Path(results_dir)
    rows = []
    for _, row in manifest.iterrows():
        result_path = results_dir / f"{row['net1']}{row['stn1']}_{row['net2']}{row['stn2']}__FastMspec.json"
        # Requiring a converged Round 1 pick to derive c_min would silently exclude exactly the
        # pairs the sweep most needs to include -- the far/currently-non-converging quartiles,
        # where convergence is ~0-4% (Round 1's own distance-quartile table). Found live while
        # testing this script: quartiles 3/4's converged pool was too small to fill 15 slots.
        # Fall back to the bare reference curve (delta=0) when no converged pick exists -- less
        # precise for that specific pair, but physically reasonable, and available for every pair
        # regardless of Round 1 outcome, which is the whole point of testing pairs that currently
        # fail.
        best_delta = 0.0
        source = "reference_curve_fallback"
        if result_path.exists():
            result = json.loads(result_path.read_text())
            if result.get("converged") and result.get("best_delta_km_s") is not None:
                best_delta = result["best_delta_km_s"]
                source = "round1_converged_pick"
        c_min = compute_c_min(best_delta, Path(ref_curve_path))
        nwh = nw_high(c_min, row["dist_km"])
        rows.append({
            "net1": row["net1"], "stn1": row["stn1"], "net2": row["net2"], "stn2": row["stn2"],
            "dist_km": row["dist_km"], "quartile": int(row["quartile"]),
            "c_min_source": source, "best_delta_km_s": best_delta, "c_min_km_s": c_min,
            "nw_high": nwh, "target_nws": [round(f * nwh, 2) for f in FRACTION_GRID],
        })
    df = pd.DataFrame(rows)

    selected = []
    for q in [1, 2, 3, 4]:
        pool = df[df["quartile"] == q]
        take = min(N_PER_QUARTILE, len(pool))
        selected.append(pool.sample(n=take, random_state=42) if len(pool) > 0 else pool)
    return pd.concat(selected, ignore_index=True)


if __name__ == "__main__":
    import sys
    manifest_csv, results_dir, ref_curve_path, out_csv = sys.argv[1:5]
    subset = select_subset(manifest_csv, results_dir, ref_curve_path)
    subset.to_csv(out_csv, index=False)
    print(f"Selected {len(subset)} pairs across {subset['quartile'].nunique()} quartiles "
          f"-> {out_csv}")
    print(subset.groupby("quartile").size())
