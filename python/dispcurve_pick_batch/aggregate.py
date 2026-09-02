"""Combines the per-work-unit JSON result files either driver produces into one manifest.csv.
Idempotent/resumable-friendly: run any time, on however many result files exist so far (a partial
batch still aggregates to a partial, valid manifest).

Usage: python3 -m dispcurve_pick_batch.aggregate <results_dir> <output_csv>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    results_dir, output_csv = Path(sys.argv[1]), Path(sys.argv[2])

    rows = []
    for f in sorted(results_dir.glob("*.json")):
        rows.append(json.loads(f.read_text()))

    if not rows:
        print(f"No result files found in {results_dir}")
        sys.exit(1)

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    print(f"Wrote {len(df)} rows to {output_csv}")
    print(f"Converged: {df['converged'].sum()}/{len(df)}")
    print(df.groupby("technique")["converged"].agg(["sum", "count"]))
    n_errors = df["error"].notna().sum()
    if n_errors:
        print(f"\n{n_errors} work units errored (not just failed to converge) -- see the "
              f"'error' column, e.g.:")
        print(df[df["error"].notna()][["net1", "stn1", "net2", "stn2", "technique", "error"]].head(10))


if __name__ == "__main__":
    main()
