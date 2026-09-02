"""Batch pipeline: recomputes all 380 Madagascar station pairs x 4 ccf_pipeline techniques
(single-taper, FastMspec, Mspec, MspecBestK), picks a phase-velocity curve for each via the
instrumented seislib picker (dispcurve_pick) scanned across the SDISPL.ASC template corridor, and
cross-validates against Sayan Swar's own precomputed MATLAB cross-spectra where available.

Two drivers exist over the same `work_unit.process()` core, to compare on real measured
throughput rather than assumption (see NOTES.md): `run_plain` (one work unit per invocation, for
a SLURM array) and `run_multiprocessing` (many work units per invocation via a process pool, for
one array task per node). `aggregate.py` combines the per-work-unit result files either produces
into one manifest.csv.
"""
