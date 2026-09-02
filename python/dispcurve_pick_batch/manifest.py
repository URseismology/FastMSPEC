"""Reads the 380-pair Madagascar station-pair catalog and expands it into the full list of
(pair, technique) work units this batch processes -- 380 x 4 = 1520.

The catalog itself (`madagascar_stn_conn_ccflist.csv`) lives on bluehive
(`PRJ_SPAC/data/test/metadata/`), not in this repo (bulk/derived data, matching this project's
established practice) -- pull via `ssh bluehive "cat <path>"` before running anything here. See
NOTES.md for the exact path and pull command.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

TECHNIQUES = ("single-taper", "FastMspec", "Mspec", "MspecBestK")


@dataclass(frozen=True)
class Pair:
    net1: str
    stn1: str
    net2: str
    stn2: str
    dist_km: float
    matched_data_path: str

    @property
    def pair_id(self) -> str:
        return f"{self.net1}{self.stn1}_{self.net2}{self.stn2}"


@dataclass(frozen=True)
class WorkUnit:
    pair: Pair
    technique: str

    @property
    def work_unit_id(self) -> str:
        return f"{self.pair.pair_id}__{self.technique}"


def load_pairs(manifest_csv: str) -> list[Pair]:
    df = pd.read_csv(manifest_csv)
    return [
        Pair(
            net1=row.net1, stn1=row.stn1, net2=row.net2, stn2=row.stn2,
            dist_km=float(row.stndist), matched_data_path=row.filelocation,
        )
        for row in df.itertuples()
    ]


def build_work_units(manifest_csv: str) -> list[WorkUnit]:
    """Technique-outer, pair-inner ordering (all single-taper work units first, then all
    FastMspec, etc.) -- NOT pair-outer -- so each technique occupies one contiguous index range.
    Mspec costs >10x every other technique (Stage 3 timing pilot: ~45 min vs. ~2.5 min cross-
    spectrum); a single SLURM array resource request sized for one technique is wrong for the
    others, so submit_plain.sbatch/submit_mp.sbatch target one technique's contiguous range at a
    time, each with its own --time/--mem sized for that technique (see TECHNIQUE_INDEX_RANGES).
    """
    pairs = load_pairs(manifest_csv)
    return [WorkUnit(pair=p, technique=t) for t in TECHNIQUES for p in pairs]


def technique_index_ranges(manifest_csv: str) -> dict[str, tuple[int, int]]:
    """{technique: (start_index, end_index)} (end exclusive) into build_work_units()'s list,
    for submitting one array job per technique with technique-appropriate resources."""
    n_pairs = len(load_pairs(manifest_csv))
    return {t: (i * n_pairs, (i + 1) * n_pairs) for i, t in enumerate(TECHNIQUES)}
