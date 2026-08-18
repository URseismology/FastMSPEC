# Octave-based verification against the original MATLAB source

[← Back to repo README](../../README.md) | See also: [python/NOTES.md](../../python/NOTES.md) (the translation this verifies) | [legacy/matlab_source/README.md](../../legacy/matlab_source/README.md)

Octave (9.4.0, `sudo apt install octave`) runs the **actual, unmodified**
`.m` files from `matlab_source/ThomsonsMethodRevisitedExperiments/` and
compares their output directly against the Python translation -- a much
stronger check than the independent-reference tests in `python/tests/`,
since this exercises the original source itself, not just a reference
reimplementation.

## What had to be supplied (not part of the original library)

- `tridisolve.m` -- the shipped file is a MEX-stub with no compiled binary
  anywhere in the codebase (see `python/NOTES.md`). Replaced with a
  hand-transcription of the Thomas algorithm from the same Golub & Van Loan
  reference its own doc comment cites -- written independently from the
  Python version, not copied from it.
- `dpss.m` (= `dpss_ref.m`) -- `octave-signal` (which would provide a real
  `dpss`) isn't installed. Supplied a from-scratch dense reference: builds
  the tridiagonal commuting matrix, diagonalizes with Octave's own `eig()`,
  and computes each eigenvalue via the textbook Toeplitz-sinc quadratic
  form `e'*B*e` -- deliberately *not* reusing `tridieig`/`tridisolve`/
  `fftfilt`, so it's an independent check, not circular.
- `datawrap.m` -- MATLAB Signal Processing Toolbox builtin, not present in
  Octave; independent shim (pad + reshape + sum), not derived from the
  Python version.
- `firstNlambdaDPSS_fixed.m` -- the original with the two bug fixes applied
  (see below), to confirm the fixes actually work in the real MATLAB
  language semantics, not just in the Python translation.

## Results (all against the same `signal.csv`, a synthetic 256-sample tone
plus noise, generated once by Python and read by both languages)

| Check | Result |
|---|---|
| `tridieig` vs Octave `eig()` on random symmetric tridiagonal | max abs diff 1.1e-15 |
| `tridisolve` vs Octave `A\b` on random symmetric tridiagonal | max abs diff 5.6e-17 |
| **Original, unmodified** `firstNlambdaDPSS.m` vs dense reference | max abs diff **0.30** -- confirms the bug is in the upstream source, not a translation artifact |
| `firstNlambdaDPSS_fixed.m` (2 fixes applied) vs dense reference | max abs diff 8.9e-16 |
| Original `transitionDPSS.m` (unmodified) vs dense reference | max abs diff 7.8e-16 |
| Original `Multitaper.m` vs Python `Multitaper` | relative L2 diff 2.1e-13 |
| Original `FastMultitaper.m` vs Python `FastMultitaper` | relative L2 diff 3.7e-16 |
| Original `MultitaperAdaptive.m` vs Python `MultitaperAdaptive` | relative L2 diff 2.1e-13, iteration count matches exactly (6) |

Conclusion: the Python translation in `python/thomson_multitaper/` matches
the actual original MATLAB source to numerical precision everywhere it was
tested, and the two bug fixes documented in `python/NOTES.md` are confirmed
against the real upstream `.m` files, not just against the translation's own
prior (buggy) behavior.

## Not covered

`transitionDPSS_modif.m` was not run here (still untested; see
`python/NOTES.md`). Only synthetic signals were used, not real seismic data.
