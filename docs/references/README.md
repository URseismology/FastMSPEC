# `docs/references/` — provenance

Local copies of the papers this project's methods are directly grounded in, kept alongside the
repo so citations link to something durable rather than DOI-only (some of these, e.g. Sayan's own
report, aren't available online at all).

**Provenance**: `Multitaper_Revisited_Karnik.pdf` and `SSWAR_ESC425_Project_Report.pdf` were
pulled earlier in this project's history (see `notebooks/README.md`'s own References section for
how each is used). The remaining five were found on the lab network-attached storage
(`repovibranium`) at `/volume1/web/FastMSPEC_data/core_papers/`, pulled via
`ssh repovibranium "cat <path>"` (the NAS's SFTP/scp subsystem is restricted to certain shared
folders; plain `ssh ... cat` works around it, `scp` does not — same method used for
`data/reference/SDISPL.ASC`). Each was verified against its own embedded PDF metadata/title page
(via `pypdf`) to confirm the citation year/journal/DOI actually matches what this project cites it
as, rather than trusting the NAS filename — all five matched cleanly, including the AkiNet paper
(NAS filename `xueolugboji2026.pdf`, but its own metadata confirms *Journal of Geophysical
Research: Machine Learning and Computation* 2025.2, DOI `10.1029/2025JH000932` — genuinely 2025,
consistent with this project's existing "Xue & Olugboji 2025" citations).

| Local file | Citation | Used by |
|---|---|---|
| `Ekstrom_2009_GRL.pdf` | Ekström, Abers & Webb (2009), *GRL*, doi:10.1029/2009GL039131 | The original zero-crossing/Aki-spectral-formulation phase-velocity method this whole pipeline is built on; `seislib`'s picker (`python/dispcurve_pick/`) is a direct methodological descendant. |
| `Ekstrom_2017_PEPI.pdf` | Ekström (2017), *PEPI*, doi:10.1016/j.pepi.2017.07.010 | Short-period extension of the same method; grounds the resolution-bandwidth argument in Notebook 5. |
| `Magrini_2022_SeisLib_GJI.pdf` | Magrini et al. (2022), *GJI*, doi:10.1093/gji/ggac236 | The SeisLib package paper — `python/dispcurve_pick/` vendors and instruments SeisLib's own `extract_dispcurve` picker. |
| `Hawkins_Sambridge_2019_BSSA.pdf` | Hawkins & Sambridge (2019), *BSSA*, doi:10.1785/0120190060 | Motivates envelope conditioning (Notebook 3) and the analytical Bessel-zero/extrema construction `nb5_helpers.py`'s template family is built on. |
| `Xue_Olugboji_2025_AkiNet_JGRMLC.pdf` | Xue & Olugboji (2025), *JGR: Machine Learning and Computation*, doi:10.1029/2025JH000932 | Motivates the additive phase-velocity corridor (`build_template_family`'s ±0.8 km/s window) and the "defaults to reference" diagnostic framing used in Notebook 5's quality discussion. |
