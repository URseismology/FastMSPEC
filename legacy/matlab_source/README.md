# Original MATLAB source

[← Back to repo README](../../README.md) | See also: [python/NOTES.md](../../python/NOTES.md) | [python/ccf_pipeline/NOTES.md](../../python/ccf_pipeline/NOTES.md) | [legacy/para_ccf_original/NOTE.md](../para_ccf_original/NOTE.md)

The MATLAB source this project's Python translation ([`../../python/`](../../python/)) was built from and verified
against, kept here so the two can be read side-by-side. Copied from the lab's Sayan Swar
Workspace codebase snapshot (dated 2026-07-13); not modified, except that only the files actually
relevant to this translation are included (the full codebase is far larger).

## Layout

- **[`ThomsonsMethodRevisitedExperiments/`](ThomsonsMethodRevisitedExperiments/)** — Santhosh Karnik's "Fast Slepian Transform" toolbox:
  DPSS/multitaper spectral estimation via a tridiagonal eigenproblem. Translated in full to
  [`python/thomson_multitaper/`](../../python/thomson_multitaper/).
- **[`lib/`](lib/)** — the ambient-noise cross-correlation (CCF) pipeline layer that calls the toolbox
  above. `ccf_compute_crosscorr_mtc_Z.m`/`_T.m` are the actual integration point with
  `FastMultitaper`; `ccf_compute_crosscorr_Z.m`/`_T.m` are the real dispatcher entry point (called
  by the only two production driver scripts found); `ccf_prepare_data_Z.m` and
  `ccf_preprocess_filter_data.m` are the SAC-loading/windowing and preprocessing stages. Translated
  to [`python/ccf_pipeline/`](../../python/ccf_pipeline/).
- **[`functions/`](functions/)** — supporting utilities the `lib/` files depend on: the `_3dim` preprocessing
  helpers, `FiltFiltM.m` (zero-phase filtering), and a minimal subset of J.M. Lilly's **jLab**
  toolbox (`jSpectral/`, `jCommon/`, `jVarfun/`) that `sleptap.m`/`mspec_fast.m` need to run.
- **[`entry_points/`](entry_points/)** — the four `a1_ccf_ambnoise_*.m` scripts originally investigated as the
  translation target, plus the one real driver script found
  (`a2_ccf_run_crosscorr_T_mdg.m`) that actually calls the `lib/` dispatcher with concrete
  parameters. Kept for context even though they weren't directly translated — see
  [`docs/plan_ccf_mtc_translation.md`](../../docs/plan_ccf_mtc_translation.md) for why the investigation moved one layer down to
  `lib/ccf_compute_crosscorr_mtc_*.m` instead.

## How to trace a translation decision back to its source

Every Python module under [`python/`](../../python/) documents, in its own `NOTES.md`
([`python/NOTES.md`](../../python/NOTES.md), [`python/ccf_pipeline/NOTES.md`](../../python/ccf_pipeline/NOTES.md)),
exactly which `.m` file(s) it corresponds to and any bugs found or judgment calls made along the
way. Start there, then open the matching file here to compare line-by-line.
