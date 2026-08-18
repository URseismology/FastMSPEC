# Plan: Translate the CCF Multitaper Cross-Correlation Pipeline (`ccf_compute_crosscorr_mtc_Z/T.m`)

[← Back to repo README](../README.md) | See also: [python/ccf_pipeline/NOTES.md](../python/ccf_pipeline/NOTES.md) (implementation status per phase) | [verification/octave_verify_ccf_pipeline/README.md](../verification/octave_verify_ccf_pipeline/README.md)

## Status update (2026-08-18, end of session)

All 6 phases below are implemented and committed. Phases 0-4 and 6's `IsMspec` dispatch path are
verified against real Octave runs of the actual, unmodified `.m` source (23/23 tests passing in
`python/ccf_pipeline/`). Found and fixed one real bug (wrong two-sided spectrum length for
certain window sizes) and found a second real MATLAB bug that crashes on every `MspecBestK` run
(a `totalMB` reference) rather than reproducing it. Phase 5 (obspy SAC loading) is implemented
and self-consistency-tested but **not yet verified against real data or Octave** -- no SAC files
were available this session. `octave-signal` got installed (thanks to you running the sudo
command remotely) enabling `ccf_butterfilt_3dim` to at least run against Octave, though its exact
match against `FiltFiltM.m`'s non-standard filtfilt variant remains an open, low-priority
discrepancy (unused in every real config found). Full details in `python/ccf_pipeline/NOTES.md`.

**What's left**: real SAC data (paths to track down are in Phase 5/6's "Required input files"
section below) to (a) verify `prepare_data.py` against Octave and real output, and (b) assemble
and verify a true end-to-end SAC-file-to-saved-output run.

## Context

The `ThomsonsMethodRevisitedExperiments` multitaper/DPSS library is already translated to Python and verified (against independent references and against the actual MATLAB source via Octave — see `~/claude-sandbox/projects/sayan-swar-translation/`). The next goal is to translate the layer of this ambient-noise cross-correlation (CCF) codebase that actually *uses* that library, so real seismic processing can move to Python.

Investigation (two rounds of code reading, summarized below) found that the four `a1_ccf_ambnoise_*.m` scripts originally named as the target do **not** call the translated library at all — they use a separate, inline `sleptap()` tapering function. The real integration point is one layer downstream: `lib/ccf_compute_crosscorr_mtc_Z.m` / `_T.m`, which is where `FastMultitaper` (from the already-translated library) is actually invoked. Per your direction, this plan targets that file, uses `obspy` in place of hand-translated SAC I/O, and accounts for `functions/jSpectral/` (the jLab-derived toolbox these files also depend on).

## What the pipeline actually does (resolved call graph)

```
ccf_prepare_data_Z.m / _T.m          (SAC loading, pairing, windowing → S1_data_mat/S2_data_mat)
        │  saves/loads via .mat, shape (day, window, samples)
        ▼
ccf_preprocess_filter_data.m         (optional: detrend, taper, prefilter, then ONE of:
                                       multitaper/FTN/OBN/spectral-whiten — mutually exclusive)
        ▼
ccf_compute_crosscorr_Z.m / _T.m     (REAL entry point — dispatcher, called by the only two
                                       actual driver scripts found: a2_ccf_run_crosscorr_T_mdg.m,
                                       ccf_synth_setup_perfrm_evl.m)
        │
        ├─ if IsMultiTaper||IsFTN||IsOBN||IsSpecWhiten: data already frequency-domain, use as-is
        ├─ elif IsMspec:  ──────────► delegates to ccf_compute_crosscorr_mtc_Z.m / _T.m  ◄── TARGET
        │                              (technique = 'Mspec' | 'FastMspec' | 'MspecBestK')
        └─ else: plain fft(data,[],3)
```

Both real example configs found in the codebase set `IsMspec=1`, so `ccf_compute_crosscorr_mtc_Z/T.m` is genuinely the live path, not a dead branch. Its three techniques:
- **`'FastMspec'`** — calls `FastMultitaper(N, Wband, cutoff, epsilon)` (your translated class) to get fused/weighted tapers, then a taper-application+FFT+cross-spectrum stage (`mspec_fast.m`'s `avgspec_xy_sayan`, a **batched, cross-spectrum generalization** of `FastMultitaper.SpectralEstimate` that doesn't exist yet in the ported library — this is the main new code needed).
- **`'Mspec'`** — classical multitaper: `sleptap()` for plain DPSS tapers, then `mspec_fast.m`'s classical path (`avgspec_xy`). `sleptap()` computes the *same* DPSS tapers your `dpss()` wrapper (backed by `scipy.signal.windows.dpss`) already produces — just via a different numerical method (MATLAB `eigs()` vs. your tridiagonal solver vs. scipy). Since only the taper *output* needs to match (not the internal method), `sleptap` does **not** need a line-by-line port — your existing `dpss()` wrapper is a legitimate substitute.
- **`'MspecBestK'`** — hybrid: `FastMultitaper` only to pick a taper count `K`, then plain `sleptap`-style tapers of that count (again → your `dpss()` wrapper).

**A resolved point of confusion worth recording**: initial analysis flagged the `N1` passed to `FastMultitaper`/`sleptap` inside `ccf_compute_crosscorr_mtc_Z.m` as possibly being a trace count rather than a sample count (which would mean the FFT axis semantics were backwards). Careful manual tracing of MATLAB's column-major `reshape(A, [], N1)` semantics indicates this is very likely **not** a bug — `N1` should still resolve to the sample count after the reshape+transpose, consistent with the real driver script's own `FastMultitaper(dt*winlength*3600+1, ...)` call (which unambiguously uses sample count). **This must be confirmed empirically, not by further reasoning** — it's exactly the kind of MATLAB indexing question that's easy to get backwards (both a prior analysis pass and a first-pass manual check on this disagreed with each other before landing here). This is Phase 0 below.

## Phased approach

### Phase 0 — Resolve the axis-semantics question empirically (Octave, synthetic data, ~30 min)
Before writing any Python: build a tiny synthetic `S1_data_mat`/`S2_data_mat` (e.g. 2 days × 3 windows × 8 samples, filled with distinguishable values), run it through the actual `ccf_compute_crosscorr_mtc_Z.m` in Octave with `technique='FastMspec'` and small dummy `Wband`/`cutoff`/`epsilon`, and print/inspect intermediate shapes (`S1_data_mat_2D` shape, what `N1` resolves to, what `mspec_fast`'s FFT operates over). This nails down the exact semantics before they get baked into a Python port. Needs the jSpectral files (`sleptap.m`, `mspec_fast.m`) and their `jCommon`/`jVarfun` utility dependencies (`vrep`, `vmean`, `fourier`, etc. — not yet fully inventoried, small follow-up read) to actually run in Octave — same "supply what's missing, run the real source" strategy that worked for the multitaper library verification.

### Phase 1 — Core engine: `ccf_compute_crosscorr_mtc_Z.m`, `'FastMspec'` technique only
The direct integration point with the already-translated library. Build:
- A Python function/class operating on synthetic `S1_data_mat`/`S2_data_mat` numpy arrays directly (no SAC/obspy yet) — mirrors the MATLAB function's reshape → `FastMultitaper` construction → taper/FFT/cross-spectrum → one-sided-to-two-sided reflection → stack.
- The new piece: a batched, cross-spectrum generalization of `FastMultitaper.spectral_estimate` (matching `avgspec_xy_sayan`'s math), added to `thomson_multitaper` or a new sibling module — reusing the existing `FastMultitaper` class for taper construction, not reimplementing DPSS math.
- Verify against Octave running the real `ccf_compute_crosscorr_mtc_Z.m` on the same synthetic input (same technique as the multitaper library's verification: independent Octave run, diff outputs), using the real example config values found (`Wband=0.001, epsilon=1e-5, cutoff=1-epsilon` from `a2_ccf_run_crosscorr_T_mdg.m`).

### Phase 2 — `'Mspec'` and `'MspecBestK'` techniques
Add the classical (non-fused) multitaper cross-spectrum path, using the existing `dpss()` wrapper in place of `sleptap()`. Verify against Octave the same way, using the second real config found (`ccf_synth_setup_perfrm_evl.m`: `technique='Mspec', Wband=9.7656e-04, cutoff=0.999, epsilon=1e-9`).

### Phase 3 — Port the `_T.m` sibling
Apply the same logic to `ccf_compute_crosscorr_mtc_T.m`, noting its two differences from Z: `'Mspec'` reads `NW_mspec`/`K_taps_mspec` directly from config instead of deriving from `Wband`; `'FastMspec'` expects a pre-built `FMTSE` object passed in (`filttype.FMTSE`) rather than constructing it inline — the Python port should support both an inline-construct and an injected-taper calling convention.

### Phase 4 — Preprocessing stage: `ccf_preprocess_filter_data.m`
Port the optional detrend/taper/prefilter/multitaper-or-FTN-or-OBN-or-whiten pipeline that produces the frequency-domain (or still-time-domain) `S1_data_mat`/`S2_data_mat` that Phase 1–3's code consumes. Includes the `_3dim` helper functions (`ccf_detrend_3dim.m`, `ccf_cos_taper_3dim.m`, `ccf_butterfilt_3dim.m`, `ccf_slepian_multitap_3dim.m`, `ccf_FTN_3dim.m`, `ccf_OBN_3dim.m`, `ccf_spectrumwhiten_smooth_3dim.m`) — all fairly mechanical numpy/scipy translations (detrend, Tukey-style taper, `scipy.signal.butter`+`filtfilt`, etc.), verified against Octave the same way.

### Phase 5 — SAC loading/windowing via obspy (replaces `ccf_prepare_data_Z.m`/`_T.m`)
Per your direction, use `obspy` (`pip install obspy`) instead of hand-translating `load_sac.m`/`readsac.m` — obspy's `read()` + `Trace.stats.sac`/`Trace.stats` map directly onto the MATLAB SAC header fields this code reads (`DELTA`, `NPTS`, `NZYEAR/NZJDAY/NZHOUR/NZMIN/NZSEC/NZMSEC` → `starttime`, `STLA/STLO/STEL`). What still needs translating (obspy doesn't give this for free): the day-file pairing-by-timestamp logic, the validation checks (matching sample rate, matching start time, minimum length, all-zero check), the min-distance station-pair filter (obspy has `gps2dist_azimuth` as a direct replacement for MATLAB's `distance()`+`deg2km()`), and the sliding-window cut-and-resample-via-interpolation logic that builds the exact `(day, window, sample)` array shape.

### Phase 6 — End-to-end wiring + real-data verification
Assemble Phases 1–5 into a Python equivalent of the `ccf_compute_crosscorr_Z.m`/`_T.m` dispatcher, run on real (or realistic synthetic) SAC data, and compare against Octave running the equivalent real `.m` files end-to-end.

## Required input files (to ask you for, or confirm synthetic data is acceptable for early phases)

Phases 0–4 need **no real data** — synthetic numpy/Octave arrays suffice and are the recommended way to validate the core math fast. Phase 5–6 need:
1. SAC waveform files, a few stations × a few days, following `{station}.{yyyy}.{jday}.{hh}.{mm}.{SS}.{COMP}.sac` — Z-component minimum; N/E components too if/when the `_T.m` path is exercised end-to-end.
2. A station list file (`sta_list.txt` / `sta_list_love.txt` format: `station lat lon elev`, whitespace-delimited).
3. An orientation CSV (`orientation.csv`: columns `net-station`, H1-orientation-degrees) — only needed for T-component rotation.
4. Confirmation of which `filttype` config to treat as "canonical" for verification — the two real configs already found in the codebase (`a2_ccf_run_crosscorr_T_mdg.m` for `FastMspec`, `ccf_synth_setup_perfrm_evl.m` for `Mspec`) are proposed as the default targets unless you have a different reference run in mind.

## Verification approach (consistent with what worked for the multitaper library)

- Synthetic-input Octave dry-runs of the actual, unmodified `.m` files at each phase, diffing against the Python port on identical input — same strategy that caught 2 real bugs and gave ~1e-12 agreement for the multitaper library.
- Any missing Octave dependencies (e.g. `jCommon`/`jVarfun` utility functions used by `mspec_fast.m`) supplied as independently-written shims where MATLAB-only, or run as-is where the real `.m` source is available (most of `jSpectral` is regular `.m` source, not MEX).
- A `tests/` suite per phase, mirroring `python/tests/test_thomson_multitaper.py`'s structure.
- `NOTES.md`-style documentation of every judgment call, flagged inconsistency (e.g. the `S2`-before-`S1` argument-order swap unique to `'FastMspec'`, which affects conjugation/sign of the cross-spectrum — must be replicated exactly), and open question, same as the existing library's `NOTES.md`.

## Not in scope for this plan (flagged, not forgotten)

- The four `a1_ccf_ambnoise_*.m` scripts themselves (SAC-to-CCF full pipeline including plotting/logging) — still not translated; may be revisited later since Phase 5's obspy loading logic overlaps significantly with what they'd need.
- `MultitaperAdaptive`-style adaptive weighting inside `mspec_fast.m` (`adaptspec`) — confirmed unused by this pipeline's actual calling pattern (always scalar `lambda`), so out of scope unless that changes.
- `sleptap_modif`, `sleptap_modif2`, `sleptap_modif3`, `mspec_loop.m`, and the three `mspec_*_bkup*.m` files — confirmed dead code (no callers in the live path) via `grep -rl`, excluded from translation scope.
- `a1_ccf_ambnoise_RTZ_NE_ADAMA.m`'s text-file "AkiEstimate" export format and its undefined `read_SACPZ`/`rm_SACPZ`/`doy2date` dependencies — out of scope until/unless that script specifically becomes a target.
