# `SDISPL.ASC` — reference Love-wave dispersion curve

Multi-mode Love-wave phase-velocity dispersion computation, Herrmann Computer-Programs-in-
Seismology ASCII format (columns: `LMODE`, `NFREQ`, `PERIOD(S)`, `FREQUENCY(Hz)`, `C(KM/S)`).

**Provenance**: found on the lab network-attached storage (`repovibranium`) at
`/volume1/web/FastMSPEC_data/madagascar_data/pre_computed_files/SDISPL.ASC`, pulled via
`ssh repovibranium "cat <path>"` (the NAS's SFTP/scp subsystem is restricted to certain shared
folders; plain `ssh ... cat` works around it, `scp` does not). This is the same file Sayan Swar's
own `phasevel_compute_slide13and14.ipynb` (also on `repovibranium`, in
`madagascar_codes/python/`) used to build the reference curve shown in his project report's
Figure 6 — filtered there, and here, to `LMODE==0` (fundamental mode) via
`[FREQUENCY_Hz, C_KM_S]`.

**Contents**: fundamental mode (`LMODE=0`) spans periods 2-2048 s (velocities 1.22-4.83 km/s,
1024 points); higher modes (`LMODE=1..8`) each start at progressively shorter periods.

**Used by**: `notebooks/05_coherence_barcode.ipynb` (`notebooks/_lib/nb5_helpers.py`), as the base
reference curve for the coherence-barcode template-matching method. Full design rationale:
`docs/coherence_barcode_design.tex` / `.pdf`.
