# Phase 0 findings: `ccf_compute_crosscorr_mtc_Z.m` empirical dry-run

[← Back to repo README](../../README.md) | See also: [python/ccf_pipeline/NOTES.md](../../python/ccf_pipeline/NOTES.md) (the translation this verifies) | [docs/plan_ccf_mtc_translation.md](../../docs/plan_ccf_mtc_translation.md) | [legacy/matlab_source/README.md](../../legacy/matlab_source/README.md)

Ran the actual, unmodified `ccf_compute_crosscorr_mtc_Z.m` (+ its real `jSpectral`/`jCommon`/
`jVarfun` dependencies, copied in here from the codebase snapshot) in Octave on a tiny synthetic
`S1_data_mat`/`S2_data_mat` (2 days x 3 windows x 8 samples, `technique='FastMspec'`) to resolve
two questions before writing any Python.

## Supplied dependencies (not part of the function under test)

Missing pieces needed to actually execute the real source in Octave, none of which change its
logic:
- `contains.m` — MATLAB string builtin not in this Octave version; minimal single-pattern shim.
- `memory_watch.m`, `allall.m`, `to_overwrite.m` — real files from the codebase (`lib/`,
  `functions/jCommon/`), just copied in; not reimplemented.
- Plus everything already verified for the multitaper library (`FastMultitaper.m`,
  `transitionDPSS.m`, `tridieig.m`, `tridisolve.m`, `dpss.m`, `datawrap.m`) and the rest of
  `jSpectral`/`jCommon`/`jVarfun` that `mspec_fast.m`/`sleptap.m` need (`mspec_fast.m`,
  `sleptap.m`, `fourier.m`, `aresame.m`, `frac.m`, `iseven.m`, `isodd.m`, `lnsd.m`, `squared.m`,
  `vfilt.m`, `vindex.m`, `vindexinto.m`, `vmean.m`, `vrep.m`, `vsize.m`, `vstd.m`, `vsum.m`).

Note: `squared.m` **does exist** in this codebase (`functions/jCommon/squared.m`) — this
retroactively resolves a caveat flagged in `python/NOTES.md` for `transitionDPSS_modif.m`, which
called `squared(e)` and was assumed to be an undefined-function bug. It isn't; the environment
used to verify the multitaper library simply didn't have this file copied in (it wasn't needed
for anything tested there). Not urgent to fix now since `transitionDPSS_modif` remains out of
scope/untested either way, but noted for accuracy.

## Finding 1: axis-semantics question — resolved, not a bug

Printed the shape of `S1_data_mat_2D = reshape(ccf_data_file.S1_data_mat, [], N1)` directly:
for the 2x3x8 test input, `size(S1_data_mat_2D) = [6, 8]` (rows = day*window = 6 traces,
columns = 8 samples), and `size(S1_data_mat_2D, 2) = 8` — the sample count, confirmed. This
matches the manual re-derivation in the plan doc and **contradicts** an earlier analysis pass
that concluded `N1` ends up being the trace count. `FastMultitaper`/`sleptap` are correctly
invoked with the window's sample length, consistent with the real driver script
(`a2_ccf_run_crosscorr_T_mdg.m`)'s own `FastMultitaper(dt*winlength*3600+1, ...)` call.
**No fix needed here** -- the Python port should use sample count for `N`, matching the MATLAB
source as literally written.

## Finding 2: real bug in the one-sided -> two-sided reflection step

`ccf_compute_crosscorr_mtc_Z.m` (and `_T.m`, same code) does:
```matlab
sig_len = length(coh_sum);      % length of the ALREADY one-sided spectrum, NOT the window length N
if mod(sig_len,2)==0
    coh_sum_neg = conj(flipud(coh_sum(2:end-1)));   % excludes DC and Nyquist
else
    coh_sum_neg = conj(flipud(coh_sum(2:end)));     % excludes only DC
end
coh_sum = [coh_sum; coh_sum_neg];
```
The correct rule for reconstructing a two-sided spectrum from a real-valued time-domain signal
depends on the parity of the *original* window length `N`, not the parity of the one-sided
spectrum's length (`floor(N/2)+1`). These two parities only agree for `N mod 4 in {1, 2}`; for
`N mod 4 in {0, 3}` the code takes the wrong branch and the reconstructed spectrum comes out
either one sample too long or one sample too short (see `check_reflection_bug.py`'s table:
confirmed for N=4,8,12 and N=7,11 -- wrong; N=6,10,10801,9,5 -- coincidentally correct).

**This is data-dependent and easy to miss**: the real production window length used in the only
concrete driver script found (`winlength=3` hr, `dt=1` -> `N=10801`) has `10801 mod 4 = 1`,
which happens to fall in the "coincidentally correct" case. The bug would bite for other window
lengths (e.g. `winlength=1` hr with `dt=1` -> `N=3601`, `3601 mod 4 = 1`, also fine; but
`winlength=2` hr -> `N=7201`, `7201 mod 4 = 1`, fine; need `dt` such that `N mod 4` lands on 0
or 3 to trigger it -- e.g. the synthetic-data example script uses `dt=0.5`, giving different N
parity behavior worth double-checking if that config is ever used for real verification).

**Plan for the Python port**: implement the reflection using the *true* window length `N`'s
parity (the numerically correct rule, equivalent to what `numpy.fft.irfft`-style reconstruction
would do), not the literal buggy MATLAB check -- flagged here and in the port's own docstring/
NOTES, same pattern as the two bugs already fixed and documented in the multitaper library.
