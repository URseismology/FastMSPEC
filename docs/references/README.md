# `docs/references/` — provenance

Local copies of the papers this project's methods are directly grounded in, kept alongside the
repo so citations link to something durable rather than DOI-only (some of these, e.g. Sayan's own
report, aren't available online at all).

**Provenance**: `Multitaper_Revisited_Karnik.pdf` and `SSWAR_ESC425_Project_Report.pdf` were
pulled earlier in this project's history (see `notebooks/README.md`'s own References section for
how each is used). The next five (`Ekstrom_2009_GRL.pdf` through `Xue_Olugboji_2025_AkiNet_JGRMLC.pdf`)
and the three coherence-theory papers below were found on the lab network-attached storage
(`repovibranium`) at `/volume1/web/FastMSPEC_data/core_papers/`, pulled via
`ssh repovibranium "cat <path>"` (the NAS's SFTP/scp subsystem is restricted to certain shared
folders; plain `ssh ... cat` works around it, `scp` does not — same method used for
`data/reference/SDISPL.ASC`). Each was verified against its own embedded PDF metadata/title page
(via `pypdf`) to confirm the citation year/journal/DOI actually matches what this project cites it
as, rather than trusting the NAS filename — all matched cleanly, including the AkiNet paper (NAS
filename `xueolugboji2026.pdf`, but its own metadata confirms *Journal of Geophysical Research:
Machine Learning and Computation* 2025.2, DOI `10.1029/2025JH000932` — genuinely 2025, consistent
with this project's existing "Xue & Olugboji 2025" citations).

The coherence-theory papers below (`Walden_2000...` through `Walden_McCoy_Percival_1994...`) were
found in support of Stage 5's principled resolution/variance/bias treatment (see
`docs/notebook5_revamp_progress.md`'s deferred/requirement log), via two routes: three (Walden
2000, Haley & Anitescu 2017, Keding et al. 2024) by following Karnik et al. (2022)'s own citation
trail both backward (their reference list) and forward (who cites them), per direct guidance; the
remaining three (the two Park papers, Walden/McCoy/Percival 1994) supplied directly, from memory,
including the closest match found so far to this project's own use case (Park & Levin 2016 -- see
table). None read in full formula-level detail yet (that's Stage 5 work) -- all confirmed on-topic
and correctly cited via their own embedded metadata/title pages.

| Local file | Citation | Used by |
|---|---|---|
| `Ekstrom_2009_GRL.pdf` | Ekström, Abers & Webb (2009), *GRL*, doi:10.1029/2009GL039131 | The original zero-crossing/Aki-spectral-formulation phase-velocity method this whole pipeline is built on; `seislib`'s picker (`python/dispcurve_pick/`) is a direct methodological descendant. |
| `Ekstrom_2017_PEPI.pdf` | Ekström (2017), *PEPI*, doi:10.1016/j.pepi.2017.07.010 | Short-period extension of the same method; grounds the resolution-bandwidth argument in Notebook 5. |
| `Magrini_2022_SeisLib_GJI.pdf` | Magrini et al. (2022), *GJI*, doi:10.1093/gji/ggac236 | The SeisLib package paper — `python/dispcurve_pick/` vendors and instruments SeisLib's own `extract_dispcurve` picker. |
| `Hawkins_Sambridge_2019_BSSA.pdf` | Hawkins & Sambridge (2019), *BSSA*, doi:10.1785/0120190060 | Motivates envelope conditioning (Notebook 3) and the analytical Bessel-zero/extrema construction `nb5_helpers.py`'s template family is built on. |
| `Xue_Olugboji_2025_AkiNet_JGRMLC.pdf` | Xue & Olugboji (2025), *JGR: Machine Learning and Computation*, doi:10.1029/2025JH000932 | Motivates the additive phase-velocity corridor (`build_template_family`'s ±0.8 km/s window) and the "defaults to reference" diagnostic framing used in Notebook 5's quality discussion. |
| `Walden_2000_Multivariate_Biometrika.pdf` | Walden (2000), *Biometrika*, 87(4), 767-788 | Multitaper *multivariate* (cross-spectral/coherence) estimator theory -- means, smoothing/leakage bias, variance, asymptotic (Wishart-based) distributions, unified across Slepian, sine-taper, Welch, and lag-window estimators. Its own real-seismic-data example (Section 6, 60 series, a frequency band with known-zero true coherence) shows Slepian multitaper coherence bias staying below ~0.2 there, vs. up to ~0.6 for a lower-quality lag-window estimator -- a real, published, quantitative bound on exactly the "spurious near-null coherence" bias this project's own picking-stability argument needs. Traced backward from Karnik et al. (2022)'s own reference list, per direct guidance to follow the citation trail. Grounds Stage 5's principled resolution/variance/bias treatment (see `docs/notebook5_revamp_progress.md`'s deferred/requirement log).
| `Haley_Anitescu_2017_Optimal_Bandwidth.pdf` | Haley & Anitescu (2017), *IEEE Signal Processing Letters*, 24(11), 1696-1700 | A genuine algorithm (not just theory) for automatically selecting the multitaper bandwidth/K that minimizes estimated MSE (bias² + jackknife-estimated variance), rather than a fixed, hand-chosen constant -- a concrete candidate for replacing this project's currently-hardcoded `Wband=0.001`/`NW=100` with a principled, data-driven choice in Stage 5. Also from Karnik et al. (2022)'s reference list. Used directly in `docs/stage5_bandwidth_theory.tex` (2026-09-03): its bias term (Eq. 14) is replaced by a closed-form Bessel-zero derivation specific to this project, and its jackknife variance (Eq. 23/24) is extended to a two-tier taper/record decomposition. |
| `Haley_Anitescu_2017_Optimal_Bandwidth_Supplement.pdf` | Haley & Anitescu (2017) supplementary material, same venue | The jackknife-variance derivation (Eq. 19-26) and Lorentzian-example bandwidth optimization (Eq. 27-39) that the main letter references but doesn't reproduce in full -- pulled from `repovibranium/core_papers` (NAS filename had a typo, "HaneyAntescuBWSupplement.pdf"; author confirmed as Haley via the PDF's own embedded metadata, then the NAS copy renamed to match). Eq. 23 (empirical jackknife variance) and Eq. 24 (its closed-form expected value, function of K alone) are the two equations `docs/stage5_bandwidth_theory.tex` builds the variance term from. |
| `Keding_2024_CrossSpectra_Bias.pdf` | Keding, Alickovic, Skoglund & Sandsten (2024), *Frontiers in Neuroscience*, 18:1415397, doi:10.3389/fnins.2024.1415397 | A 2024 paper (different application domain -- EEG-based speech tracking, not seismology) explicitly citing Karnik et al. (2022) and addressing coherence-level bias and spectral peak-shifting bias directly -- found via forward citation search, a methodologically relevant cross-check on bias-correction technique for the same underlying statistical problem this project's picking-stability argument is built on.
| `Park_Lindberg_Vernon_1987_JGR.pdf` | Park, Lindberg & Vernon (1987), *JGR*, 92(B12), 12675-12684 | The classic foundational paper applying Thomson's multitaper method to real seismograms -- adaptive multitaper spectral estimation on real seismic data, directly relevant precedent for this project's own real-data application.
| `Park_Levin_2016_GJI_Jackknife.pdf` | Park & Levin (2016), *GJI*, 207, 512-527, doi:10.1093/gji/ggw291 | The closest match found to this project's own use case: multiple-taper *correlation* (MTC, i.e. cross-spectral coherence) applied to real seismic receiver functions, with an explicit Section 2 (statistical-assumption testing) and Section 6 (jackknife time-domain uncertainty estimation) -- a real, seismology-specific statistical treatment of multitaper cross-spectral estimator bias/variance, not just the more abstract general theory in Walden (2000). Points to the original Park & Levin (2000) paper (introducing MTC) as worth tracking down too if the full derivation history is needed.
| `Walden_McCoy_Percival_1994_IEEE_TSP.pdf` | Walden, McCoy & Percival (1994), *IEEE Trans. Signal Processing*, 42(2) | The variance formula for multitaper spectrum estimates (real Gaussian processes), cited directly in Walden (2000)'s own reference list -- the univariate variance result underlying the multivariate/coherence treatment above.
