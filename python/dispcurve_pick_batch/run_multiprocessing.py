"""Multiprocessing driver: processes a SLICE of work units (many, not one) per invocation, fanned
out across a `multiprocessing.Pool` of worker processes -- the in-node-parallelism counterpart to
`run_plain.py`'s one-work-unit-per-invocation model. Matches one SLURM array task claiming a full
node and using all its cores, instead of 1520 separate array tasks each paying fresh
Python/numpy/scipy/obspy import overhead for a single work unit. See NOTES.md "Plain vs.
multiprocessing" for why this pair of drivers exists and how to compare them.

Usage (run as a module, from the `python/` directory):
    python3 -m dispcurve_pick_batch.run_multiprocessing <manifest_csv> \
        <results_dir> <start_index> <end_index> [--workers N]

(No ref_curve_path argument since Stage 4.5: work_unit.process() now builds each pair's own hybrid
ADAMA+GDM52 reference curve internally from that pair's coordinates -- see work_unit.py's module
docstring.)

Processes work units [start_index, end_index) (Python-style half-open range), skipping any whose
result file already exists (same idempotent-resume behavior as run_plain.py).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import os
from pathlib import Path

from .manifest import build_work_units
from .work_unit import process


def _process_one(args):
    wu, results_dir = args
    out_path = results_dir / f"{wu.work_unit_id}.json"
    if out_path.exists():
        return wu.work_unit_id, "skipped"
    result = process(wu.pair, wu.technique)
    out_path.write_text(json.dumps(result.as_dict(), indent=2))
    return wu.work_unit_id, f"converged={result.converged} runtime={result.runtime_s:.1f}s error={result.error}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest_csv")
    parser.add_argument("results_dir")
    parser.add_argument("start_index", type=int)
    parser.add_argument("end_index", type=int)
    parser.add_argument("--workers", type=int, default=None,
                         help="default: os.cpu_count() (SLURM sets this correctly via cgroup limits when "
                              "the task requests --cpus-per-task, so the plain os.cpu_count() default is "
                              "usually right without extra SLURM-env-var plumbing)")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    work_units = build_work_units(args.manifest_csv)
    sliced = work_units[args.start_index:args.end_index]
    n_workers = args.workers or os.cpu_count()
    print(f"Processing {len(sliced)} work units [{args.start_index}, {args.end_index}) "
          f"with {n_workers} worker processes")

    tasks = [(wu, results_dir) for wu in sliced]
    with multiprocessing.Pool(n_workers) as pool:
        for work_unit_id, status in pool.imap_unordered(_process_one, tasks):
            print(f"{work_unit_id}: {status}")


if __name__ == "__main__":
    main()
