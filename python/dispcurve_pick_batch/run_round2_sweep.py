"""Round 2 driver: runs FastMspec at several candidate bandwidths per pair, for the stratified
NW-sweep subset selected by round2_prep.py. One work unit per SLURM array task (matches
run_plain.py's model -- no shared-node contention, each sweep point gets its own dedicated time
budget). See docs/notebook5_revamp_progress.md's "Round 2 design" section for the full rationale.

Usage (run as a module, from the `python/` directory):
    python3 -m dispcurve_pick_batch.run_round2_sweep <subset_csv> <main_manifest_csv> \
        <ref_curve_path> <results_dir> <work_unit_index>

<subset_csv>: round2_prep.py's output (60 pairs, each carrying its own nw_high). Target NW values
are recomputed directly from nw_high * FRACTION_GRID here, not parsed from the subset CSV's own
`target_nws` column (that column round-trips through pandas' str() of a Python list -- fragile to
re-parse; nw_high is a clean float column and FRACTION_GRID is a shared constant, so recomputing
is both simpler and exact).
<main_manifest_csv>: the standard 380-pair catalog, needed only for matched_data_path (round2_prep
's own subset CSV doesn't carry it).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from .manifest import Pair
from .round2_prep import FRACTION_GRID, N_SAMPLES
from .work_unit import process


def build_sweep_items(subset_csv: str, main_manifest_csv: str) -> list[tuple[Pair, float]]:
    """Returns [(Pair, target_nw), ...], subset row-major then fraction-inner -- 60 pairs x 5
    fractions = 300 items, in the same deterministic order every time (pandas preserves row
    order; FRACTION_GRID is a fixed module-level list).
    """
    subset = pd.read_csv(subset_csv)
    main = pd.read_csv(main_manifest_csv).rename(columns={"stndist": "dist_km"})
    lookup = {
        (r.net1, r.stn1, r.net2, r.stn2): r.filelocation
        for r in main.itertuples()
    }

    items = []
    for row in subset.itertuples():
        key = (row.net1, row.stn1, row.net2, row.stn2)
        matched_data_path = lookup.get(key)
        if matched_data_path is None:
            raise KeyError(f"Pair {key} from subset_csv not found in main_manifest_csv")
        pair = Pair(net1=row.net1, stn1=row.stn1, net2=row.net2, stn2=row.stn2,
                    dist_km=float(row.dist_km), matched_data_path=matched_data_path)
        for frac in FRACTION_GRID:
            target_nw = round(frac * row.nw_high, 2)
            items.append((pair, target_nw))
    return items


def main():
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)
    subset_csv, main_manifest_csv, ref_curve_path, results_dir, index_str = sys.argv[1:6]
    index = int(index_str)

    items = build_sweep_items(subset_csv, main_manifest_csv)
    if not (0 <= index < len(items)):
        raise IndexError(f"work unit index {index} out of range [0, {len(items)})")
    pair, target_nw = items[index]

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"{pair.pair_id}__FastMspec__NW{target_nw:.2f}.json"
    if out_path.exists():
        print(f"skip (already done): {out_path.name}")
        return

    wband = target_nw / N_SAMPLES
    result = process(pair, "FastMspec", Path(ref_curve_path), wband_override=wband)
    out_path.write_text(json.dumps(result.as_dict(), indent=2))
    print(f"{out_path.name}: NW={target_nw} wband={wband:.6f} converged={result.converged} "
          f"runtime={result.runtime_s:.1f}s error={result.error}")


if __name__ == "__main__":
    main()
