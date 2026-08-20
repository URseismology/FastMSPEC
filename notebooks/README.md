# Documentation notebooks

[← Back to repo README](../README.md)

Four pre-run Jupyter notebooks connecting the theory behind this repo's `thomson_multitaper` and
`ccf_pipeline` packages to why they matter for Sayan's ambient-noise cross-correlation problem, and
to real data. Each is committed with its outputs already saved, so it renders fully on GitHub
without needing to be re-executed.

| Notebook | What it covers | Status |
|---|---|---|
| [`01_multitaper_theory.ipynb`](01_multitaper_theory.ipynb) | Why multitaper spectral estimation, why "Fast" — reproduces Figs 1-8 and Tables I-II of Karnik, Romberg & Davenport (2022) | Figs 1-3/Table I are fresh implementations (no MATLAB driver exists for them anywhere in this codebase); Figs 4-8/rest of Table II are ported from `Comparison_ARMA_largescale_*.m`/`FastMultitaper_SpeedTest.m`, **at reduced scale** — see the notebook's own scale-note callouts for why (the paper's exact N=2^18/1000-trial parameters were benchmarked and found computationally infeasible here) |
| [`02_why_cross_spectra.ipynb`](02_why_cross_spectra.ipynb) | Why this matters for ambient-noise cross-correlation specifically: the Bessel-coherence link to phase velocity, the Love-wave problem, the scale problem, and the bridge from Notebook 1's auto-spectra to the pipeline's cross-spectra | New synthesis, grounded in Sayan's report/presentation and this repo's own `fast_cross_spectrum.py` |
| [`03_fastmspec_application.ipynb`](03_fastmspec_application.ipynb) | The pipeline run on real data: SA53/SA58 (all three `IsMspec` techniques compared), MTAN/RUNG (Love-wave SNR), a synthetic NLNM stability demo, and real dispersion-curve picking via `seislib` | New real-data results — **includes two honest negative/anomalous findings** (see below), not smoothed over |
| [`04_coda_correlation_future_work.ipynb`](04_coda_correlation_future_work.ipynb) | Roadmap for coda-correlation, explicitly out of scope for Sayan's course project | Scaffold only, no MATLAB or Python implementation exists to port |

## Worth knowing before reading Notebook 3

Notebook 3 documents an honest, non-trivial investigation, not a straight-line confirmation of the
report's numbers. Summary, fullest detail in the notebook itself (Sections 2, 2b, and 4):

- **MTAN/RUNG's SNR anomaly is real and pair-specific, not a rotation-code bug.** Section 2 first
  found single-taper SNR *beating* FastMspec on MTAN/RUNG's transverse (Love-wave) component — the
  opposite of the report's numbers. Section 2b widened the check to a second station pair
  (SA58/SA53) and both components (transverse and vertical): **SA58/SA53 shows FastMspec winning
  clearly on both components**, with numbers close to the report's own (report: 15.2→18.3 and
  5.7→10.3 dB; SA58/SA53 here: 7.4→15.3 dB (Z) and 12.4→19.5 dB (T)) — while **MTAN/RUNG shows the
  same anomaly on its vertical component too**, which needs no rotation at all. This rules out the
  new N/E→R/T rotation code as the cause and points instead at something specific to the MTAN/RUNG
  dataset (the network's earlier 1994 "XD" deployment, vs. SA53/SA58's later 1998 "XA" one).
- **A frequency-domain, Bessel-model-based SNR alternative** (`bessel_fit_quality`, grid-searching a
  best-fit phase velocity against Aki's $J_0$ coherence prediction) was added per your suggestion —
  it doesn't yet cleanly separate the two methods either (raw RMS residual is dominated by overall
  coherence strength, not Bessel-shape quality; needs amplitude normalization to be a fair metric).
- **Envelope conditioning** (motivated by, not literally ported from, Hawkins & Sambridge 2019 and
  Xue & Olugboji 2025's AkiNet — see References) surfaced a genuinely useful finding: MTAN/RUNG's
  NCF envelope peaks at 34 s lag, implying a physically plausible ≈3.2 km/s arrival, not a
  near-zero-lag artifact as first hypothesized — a real improvement over `calc_snr_onesided`'s fixed
  window — but it still didn't resolve the anomaly or favor FastMspec in the Bessel-fit metric.
- **Dispersion-curve picking** (Section 4) did not converge for *either* method on MTAN/RUNG
  (max |coherence| ≈ 0.21), a different outcome than the report's Fig. 6-7 — plausibly the same
  underlying low-coherence issue as the SNR anomaly, on this same pair.

Section 1 (SA53/SA58, `IsMspec` techniques on the already end-to-end-verified pipeline path) is not
affected by any of this. The open item worth raising with Sayan: what's different about the
MTAN/RUNG dataset or his exact original processing that this notebook's fresh rotation/SNR/Bessel
code isn't capturing.

## Running these yourself

```bash
cd python && python3 -m pip install -r requirements.txt
python3 -m pip install jupyter nbformat nbclient matplotlib pandas seislib
cd ../notebooks
jupyter nbconvert --to notebook --execute --inplace 01_multitaper_theory.ipynb  # ~15-20 min
jupyter nbconvert --to notebook --execute --inplace 02_why_cross_spectra.ipynb  # ~1 min
jupyter nbconvert --to notebook --execute --inplace 03_fastmspec_application.ipynb  # ~5-10 min, needs data/ (see repo root README)
```

`seislib` (Notebook 3, Section 4) needs `python3-dev` (or your platform's equivalent Python
development headers) installed system-wide to compile a Cython extension during `pip install` —
if that install fails with a `Python.h: No such file or directory` error, that's the fix.

`_lib/` holds the build scripts (`build_nb1.py`..`build_nb4.py`) that generate these notebooks from
scratch via `nbformat`, plus shared helper code (`karnik_figures.py`, `nb3_helpers.py`) — useful if
you want to see exactly how a figure was produced, or to regenerate a notebook after editing its
build script.

## References

Full source documents are copied into [`../docs/references/`](../docs/references/) for direct
linking (Sayan's report and presentation are unpublished coursework, not available online, so
that's the only way to link them; the two published papers are linked to their DOIs below as the
primary source, with the local copy as a convenience mirror).

- Karnik, S., Romberg, J., & Davenport, M. A. (2022). Thomson's multitaper method revisited.
  *IEEE Transactions on Information Theory*, 68(7), 4864-4891.
  [doi:10.1109/TIT.2022.3151415](https://doi.org/10.1109/TIT.2022.3151415) —
  [local copy](../docs/references/Multitaper_Revisited_Karnik.pdf). The paper Notebook 1 reproduces
  figures from.
- Karnik, S., Zhu, Z., Wakin, M. B., Romberg, J., & Davenport, M. A. (2019). The fast Slepian
  transform. *Applied and Computational Harmonic Analysis*, 46(3), 624-652.
  [doi:10.1016/j.acha.2017.07.005](https://doi.org/10.1016/j.acha.2017.07.005). The earlier paper
  this repo's whole `thomson_multitaper` library and name are built on.
- Thomson, D. J. (1982). Spectrum estimation and harmonic analysis. *Proceedings of the IEEE*,
  70(9), 1055-1096.
  [doi:10.1109/PROC.1982.12433](https://doi.org/10.1109/PROC.1982.12433). The original multitaper
  method both papers above revisit/accelerate.
- Swar, S. K. (2025). *Improving Ambient Noise Correlations with Multi-Taper Method* (ESC 425
  Planetary Seismoacoustics course project report, advised by T. Olugboji). University of
  Rochester. Unpublished. [local copy](../docs/references/SSWAR_ESC425_Project_Report.pdf) |
  [presentation slides](../docs/references/SSWAR_ESC425_Presentstion.pptx). The report Notebooks
  2-3 are directly grounded in and validate.
- Olugboji, T., & Xue, S. (2022). A short-period surface-wave dispersion dataset for model
  assessment of Africa's crust: ADAMA. *Seismological Research Letters*, 93(3), 1943-1959.
  [doi:10.1785/0220210355](https://doi.org/10.1785/0220210355). Source of the Madagascar station
  network and the Rayleigh-vs-Love Bessel-fit observation Notebook 2 builds on.
- Aki, K. (1957). Space and time spectra of stationary stochastic waves, with special reference to
  microtremors. *Bulletin of the Earthquake Research Institute*, 35, 415-456. (No DOI; pre-dates
  DOI registration.) Origin of the coherence-to-Bessel-function relationship Notebook 2 Section 1
  and Notebook 3's Bessel-fit diagnostic are built on.
- Yokoi, T., & Margaryan, S. (2008). Consistency of the spatial autocorrelation method with
  seismic interferometry and its consequence. *Geophysical Prospecting*, 56(3), 435-451.
- Ekström, G., Abers, G. A., & Webb, S. C. (2009). Determination of surface-wave phase velocities
  across USArray from noise and Aki's spectral formulation. *Geophysical Research Letters*, 36(18).
  [doi:10.1029/2009GL039131](https://doi.org/10.1029/2009GL039131).
- Ekström, G. (2014). Love and Rayleigh phase-velocity maps, 5-40 s, of the western and central USA
  from USArray data. *Earth and Planetary Science Letters*, 402, 42-49.
  [doi:10.1016/j.epsl.2013.11.022](https://doi.org/10.1016/j.epsl.2013.11.022).
- Ekström, G. (2017). Short-period surface-wave phase velocities across the conterminous United
  States. *Physics of the Earth and Planetary Interiors*, 270, 168-175.
  [doi:10.1016/j.pepi.2017.07.010](https://doi.org/10.1016/j.pepi.2017.07.010). These three,
  together with Aki (1957) above, are the literature precedent Notebook 3's Section 1f
  bandwidth/zero-crossing-spacing argument is grounded in -- large-scale, noise-derived
  phase-velocity work built on the same Aki spectral formulation.
- Tkalčić, H., Phạm, T. S., & Wang, S. (2020). The Earth's coda correlation wavefield: Rise of the
  new paradigm and recent advances. *Earth-Science Reviews*, 208, 103285.
  [doi:10.1016/j.earscirev.2020.103285](https://doi.org/10.1016/j.earscirev.2020.103285). Notebook
  4's entire theoretical grounding.
- Magrini, F., Lauro, S., Kästle, E., & Boschi, L. (2022). Surface-wave tomography using SeisLib:
  a Python package for multiscale seismic imaging. *Geophysical Journal International*, 231(2),
  1011-1030. [doi:10.1093/gji/ggac236](https://doi.org/10.1093/gji/ggac236). The `seislib` package
  Notebook 3 Section 4 uses (referred to as "SeisLab" in Sayan's report — the same package,
  confirmed by cross-checking author/citation).
- Peterson, J. R. (1993). *Observations and modeling of seismic background noise* (Open-File
  Report 93-322). U.S. Geological Survey. Source of the New Low Noise Model Notebook 3 Section 3's
  synthetic test is built on (via `obspy`'s built-in implementation, not hand-transcribed).
- Hawkins, R., & Sambridge, M. (2019). An adjoint technique for estimation of interstation phase
  and group dispersion from ambient noise cross correlations. *Bulletin of the Seismological
  Society of America*, 109(5), 1716-1728.
  [doi:10.1785/0120190060](https://doi.org/10.1785/0120190060). Motivates (not literally
  implemented — the full paper was paywalled from this environment) Notebook 3's envelope
  conditioning for the Bessel-fit dispersion diagnostic.
- Xue, S., & Olugboji, T. (2025). AkiNet: A physics-informed AI for wave extraction from noise.
  *Journal of Geophysical Research: Machine Learning and Computation*, 2(4), e2025JH000932.
  [doi:10.1029/2025JH000932](https://doi.org/10.1029/2025JH000932). Same motivating role as
  Hawkins & Sambridge above, for the same reason (paywalled).
