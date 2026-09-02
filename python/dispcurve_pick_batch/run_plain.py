"""Plain driver: processes exactly ONE work unit per invocation, selected by index. Matches a
SLURM array task model 1:1 (`--array=0-1519`, `SLURM_ARRAY_TASK_ID` as the index) -- the
no-multiprocessing baseline `submit_plain.sbatch` drives, to compare against
`run_multiprocessing.py`'s in-node worker-pool approach on real measured throughput (see NOTES.md
"Plain vs. multiprocessing" for the comparison this pair of drivers exists to run).

Each work unit's result is written as its own small JSON file (not appended to a shared CSV --
1520 concurrent SLURM tasks writing to one file is a race condition waiting to happen); run
`aggregate.py` afterward to build the final manifest.csv from all of them.

Usage (run as a module, from the `python/` directory, so the relative imports below resolve --
NOT `python3 run_plain.py`):
    python3 -m dispcurve_pick_batch.run_plain <manifest_csv> <ref_curve_path> <results_dir> <work_unit_index>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from .manifest import build_work_units
from .work_unit import process


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    manifest_csv, ref_curve_path, results_dir, index_str = sys.argv[1:5]
    index = int(index_str)

    work_units = build_work_units(manifest_csv)
    if not (0 <= index < len(work_units)):
        raise IndexError(f"work unit index {index} out of range [0, {len(work_units)})")
    wu = work_units[index]

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"{wu.work_unit_id}.json"
    if out_path.exists():
        print(f"skip (already done): {wu.work_unit_id}")
        return

    result = process(wu.pair, wu.technique, Path(ref_curve_path))
    out_path.write_text(json.dumps(result.as_dict(), indent=2))
    print(f"{wu.work_unit_id}: converged={result.converged} runtime={result.runtime_s:.1f}s "
          f"error={result.error}")


if __name__ == "__main__":
    main()
