# Translation notes: ThomsonsMethodRevisitedExperiments -> Python

[← Back to repo README](../README.md) | See also: [verification/octave_verify_multitaper/README.md](../verification/octave_verify_multitaper/README.md) (how this was verified against real MATLAB) | [legacy/matlab_source/ThomsonsMethodRevisitedExperiments/](../legacy/matlab_source/ThomsonsMethodRevisitedExperiments/) (the original .m files)

Source: `PRJ_SPAC/codes/test/matlab/functions/ThomsonsMethodRevisitedExperiments/`
(the live/current directory in the downloaded codebase snapshot, not the
`bkup/codes_020212025` or `bkup/codes_02202025` copies -- see "Which source
directory" below). Called by `a1_ccf_ambnoise_*.m`, `functions/jSpectral/*.m`,
and the `synthetics/` plotting scripts, per the ambient-noise cross-correlation
pipeline in this codebase.

Status: **translated and numerically verified against independent references**
(not against the original MATLAB directly -- no MATLAB/Octave was available on
this machine). `python3-venv`/`pip` weren't installable via `apt` without a
sudo password, so `pip` was bootstrapped user-locally via
`https://bootstrap.pypa.io/get-pip.py --user --break-system-packages`, then
`numpy`/`scipy` installed the same way (see `requirements.txt`).

`tests/test_thomson_multitaper.py` (9 tests, all passing) checks:
- `tridieig`/`tridisolve` against `numpy.linalg.eigvalsh` / `numpy.linalg.solve`
  on small dense symmetric (tridiagonal) systems -- exact to ~1e-10.
- `transition_dpss` eigenvalues against `scipy.signal.windows.dpss` across
  three (N, W) configurations -- exact to ~1e-6.
- `first_n_lambda_dpss` against `scipy.signal.windows.dpss` across three
  configurations -- exact to ~1e-15 **after fixing two bugs found in the
  original `firstNlambdaDPSS.m`** (below) -- before the fixes, errors were as
  large as 0.15 in the eigenvalue.
- `Multitaper`, `MultitaperAdaptive`, and `FastMultitaper` all correctly
  recover a known tone frequency from a synthetic noisy sinusoid, agree with
  each other to ~1e-15 relative error, and run cleanly in both the M=N and
  M<N (datawrap/decimated) code paths and both cutoff modes (integer K vs.
  fractional eigenvalue threshold).

`transition_dpss_modif` remains untested -- see its section below.

Run the suite yourself: `cd python && PYTHONPATH=. python3 tests/test_thomson_multitaper.py`
(or `python3 -m pytest tests/ -v` if pytest is installed).

## Which source directory

The tarball contains four copies of this library:
- `PRJ_SPAC/codes/test/matlab/functions/ThomsonsMethodRevisitedExperiments/` --
  the one used here. It's the only copy with `firstNlambdaDPSS.m` and
  `transitionDPSS_modif.m`, i.e. it looks like the most current version, and
  it's under the live `codes/` tree (not a `bkup/`), matching the wiki's
  directory listing and the `a1_ccf_ambnoise_*.m` callers.
- `PRJ_SPAC/bkup/codes_020212025/...` -- an exact name match for what was
  asked for ("codes_020212025"), but it's a dated backup snapshot missing the
  two newer files above.
- `PRJ_SPAC/bkup/codes_02202025/...` -- an older backup (Feb 20 2025).
- `PRJ_SPAC/reference/thompson_revisit_karnik/...` -- looks like the original
  upstream Karnik toolbox with extra ARMA-comparison scripts not used
  elsewhere in this codebase.

If the live `codes/` version isn't actually the one you want translated,
say so and I'll redo it against `bkup/codes_020212025` instead -- the core
functions are identical except for the two additions.

## Per-file notes

**tridisolve.m** -- ships as a MATLAB MEX-file stub with no visible source
(just a doc comment citing Golub & Van Loan, "Matrix Computations" 2nd ed.,
p.156). Implemented as the standard Thomas algorithm for a symmetric
tridiagonal system, per that reference. Not literally transcribed (there was
nothing to transcribe) -- verify against a known tridiagonal solve if
possible.

**tridieig.m** -- full source was available; translated near-literally,
including keeping 1-based-style index bookkeeping via padded arrays, to
minimize risk in this Sturm-sequence bisection algorithm. This one is worth
a direct numerical spot-check against `numpy.linalg.eigvalsh` on a small
tridiagonal test matrix.

**transitionDPSS.m** -- translated with MATLAB loop/index *values* preserved
exactly (1-based conceptually), converting to 0-based only at the point of
array access. `fftfilt(flipud(e), e)' * s` (used to get each taper's
concentration eigenvalue without a full quadratic form) is reproduced as
`np.convolve(b, x, mode='full')[:n]` dotted with the sinc kernel -- this
matches MATLAB's `fftfilt` semantics when both signals have equal length
`n`, but hasn't been checked against a real MATLAB run.

One structural note carried over as-is, not a translation choice: after the
inner eigenvector loop, the source unconditionally overwrites `index0` (or
`index1`) if the search range reached the end of the valid range (`k2==n` /
`k1==1`) -- even if the loop had already `break`-ed early on the epsilon
threshold. In that edge case the .m source's own index bookkeeping ends up
inconsistent (`index0=idx+j` vs. `index0=idx+j-1` from the break branch).
Preserved verbatim rather than "fixed" since I don't know which behavior
downstream code expects; flagging in case it matters.

**transitionDPSS_modif.m** -- lowest confidence translation. The MATLAB
source has two apparent bugs, both addressed here rather than reproduced,
since reproducing them would just mean re-implementing a MATLAB error:
- `squared(e)` is not a MATLAB builtin (line ~241 in the .m). Translated as
  `e.^2` / the L2 norm, matching what a "sum of squares" normalization
  clearly intends.
- `S1b` is preallocated as `zeros(Nold,0)` even in the branch that goes on
  to write `k2-k1+1` columns into it -- would error in MATLAB. Preallocated
  with the correct column count instead, mirroring `S0b`'s equivalent branch.
This function also isn't called by anything else in this directory or by the
`a1_ccf_ambnoise_*.m` callers found so far -- worth confirming it's actually
needed before investing more verification effort in it.

**firstNlambdaDPSS.m** -- the `if(k(1) > 1)` branch that would extend the
search range is unreachable given the fixed `k=[1,firstN]` starting range
(dead code in the original), so it's omitted rather than translated. Two
real bugs found and fixed here, both confirmed by comparing against
`scipy.signal.windows.dpss` (see Status above):
1. Missing `v = v(end:-1:1);` before using `v(j)` as the inverse-iteration
   shift (present in the analogous block of `transitionDPSS.m`, absent
   here) -- without it, each shift is paired with the wrong tridiagonal
   eigenvalue.
2. The symmetric/antisymmetric polarization check uses
   `mod(k(2)+1-j,2)` (the *descending*-index convention from
   `transitionDPSS.m`'s upper-region loop) while the seed vector uses the
   *ascending* `(j+k(1)-1)*t` (from that file's lower-region loop) -- these
   only agree by coincidence, and mismatched they select the wrong
   symmetrize branch, which can badly corrupt an already-converged
   eigenvector. Fixed by using the same index for both.

**Multitaper.m** -- straightforward, uses `scipy.signal.windows.dpss` as the
drop-in for MATLAB's `dpss()`. `datawrap` (time-domain aliasing before an
M-point FFT when the taper length N exceeds M) is implemented via
pad-then-reshape-then-sum.

**MultitaperAdaptive.m** -- two deviations from the source, both noted
inline in the module docstring:
1. The `cutoff > eps` branch calls `median(1, ceil(...)+10, N)` without the
   `[...]` that `Multitaper.m` uses for the equivalent line -- almost
   certainly a typo, since 3-arg `median()` doesn't mean "median of these
   three numbers" in MATLAB. Translated as the evident intent.
2. Unlike `Multitaper.m`, this file never clamps `K = max(K, 1)`. Added the
   clamp anyway since an unclamped `K=0` would crash the estimator.

**FastMultitaper.m** -- translated directly on top of `transition_dpss`.

## Octave verification (2026-08-18)

Octave 9.4.0 was installed and used to run the **actual original .m files**
directly, comparing against the Python port on identical input. Full
methodology and results table in `../verification/octave_verify_multitaper/README.md`. Summary:
every function tested (tridieig, tridisolve, transitionDPSS, the fixed
firstNlambdaDPSS, Multitaper, FastMultitaper, MultitaperAdaptive) matches
the Python translation to ~1e-12 or better, and running the **unmodified**
`firstNlambdaDPSS.m` in Octave reproduced the same 0.30 max-error pattern
found in the initial (pre-fix) Python translation -- confirming both bugs
are real defects in the upstream MATLAB source, not translation mistakes.
This resolves the "unverified against MATLAB directly" caveat from the
previous status update, for everything except `transitionDPSS_modif.m` and
real (non-synthetic) data.

## MEX dependency audit

Searched the full codebase tarball for `.mex*`/`.dll`/`.so`/`.c`/`.cpp`/`.oct`
files. The only one relevant to this library is `tridisolve.m` itself, whose
"implementation" is a MEX-file stub with no compiled binary shipped at all (no
`.mexw64` anywhere in the archive) -- just the doc comment. Already handled:
translated from the algorithm the docstring itself cites, and verified to
~1e-10 against `numpy.linalg.solve`. No other file in
`ThomsonsMethodRevisitedExperiments/` has a MEX or native-code dependency.
(A `FilterX.c` exists elsewhere in `functions/`, but it backs a sibling
function outside this library, not something this translation depends on.)

## Getting MATLAB (or Octave) for direct verification

Real MATLAB is commercial/license-gated and can't be installed the way
everything else here has been -- it needs you to authenticate to MathWorks
with a licensed account (UR almost certainly has a campus license) and
download/activate it yourself. GNU Octave is a free, no-license alternative
that can run these `.m` files directly for an output diff against the Python
translation -- the strongest remaining verification step, since everything
so far has only been checked against independent Python references
(scipy/numpy), not the original MATLAB. Installing Octave needs
`sudo apt install octave` (the one step needing an interactive password).

## What's genuinely unverified

- `transition_dpss_modif` (the N=512-reference / spline-interpolation
  variant) -- not exercised by the test suite. It also isn't called by
  anything else in this directory or by the `a1_ccf_ambnoise_*.m` callers
  found so far, so it's untested pending confirmation it's actually needed.
- Everything here has been validated against *independent* Python
  references (dense eigensolvers, scipy's DPSS), not against the original
  MATLAB on identical input. If MATLAB/Octave ever becomes available on
  this machine, running the original `.m` files on the same synthetic
  signal and diffing outputs numerically would be the strongest remaining
  check, particularly for `MultitaperAdaptive`'s iterative weighting (only
  checked indirectly via its peak-frequency output, not its exact
  convergence trajectory).
- Real seismic data hasn't been run through this yet -- only synthetic
  tones. Worth a pass on an actual ambient-noise cross-correlation segment
  from this project before treating it as a drop-in replacement for the
  `.m` pipeline.
