# Jackknife variance closed-form check

[← Back to repo README](../../README.md) | See also: [docs/stage5_bandwidth_theory.tex](../../docs/stage5_bandwidth_theory.tex) | [docs/notebook5_revamp_progress.md](../../docs/notebook5_revamp_progress.md) (Stage 5 log, 2026-09-03)

Checks whether Haley & Anitescu (2017)'s closed-form expected jackknife variance (supplement
Eq. 24, a function of the taper count `K` alone via the trigamma function) is a safe proxy for
the literal delete-one jackknife (Eq. 23) when applied to a **cross-spectrum** rather than the
auto-spectrum it was derived for -- the variance term `docs/stage5_bandwidth_theory.tex` proposes
using as FastMspec's cheap, non-per-taper-materializing variance estimate.

## Getting the data (not committed -- matches this repo's established practice)

Same file used by `verification/skrh_band_real_data/`:
```bash
ssh bluehive "cat /scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/test/processed_data/love/madagascar/AFSKRH_XVBAND_win_3_all_matched_data.mat" > AFSKRH_XVBAND_win_3_all_matched_data.mat
```
Then: `python3 validate_jackknife_closed_form.py` (or pass the `.mat` path as an argument).

## What it checks, and what it found

For `K = 5, 13, 33`, computes the real per-taper cross-spectra (kept individually, not summed --
the thing FastMspec's architecture deliberately avoids materializing) via Mspec's classical DPSS
approach, then compares the literal Eq. 23 jackknife variance against the Eq. 24 closed-form
prediction for the same `K`.

**They don't agree, and the disagreement grows with `K` rather than shrinking**: the
literal-vs-closed-form ratio is 1.9x at `K=5`, 1.5x at `K=13`, 8.3x at `K=33`. Binning the
literal jackknife variance by per-taper magnitude (a proxy for proximity to a coherence null)
at `K=33` shows why: the lowest-magnitude quartile has jackknife variance roughly 20-30x the
middle quartiles (0.037 vs. 0.001-0.002).

**Root cause**: Eq. 24's derivation assumes each per-taper quantity is `chi^2(2)`-distributed --
true for an auto-spectrum (the squared magnitude of one complex Gaussian) but not for a
cross-spectrum term (the product of *two* correlated complex Gaussians), whose distribution near
a coherence null is qualitatively different -- the same `~1/K` near-null instability documented
in Walden (2000). Eq. 24 systematically underestimates variance exactly where a picker's
crossing decisions are made.

**Consequence for Stage 5** (see `docs/stage5_bandwidth_theory.tex` Section 5.1's 2026-09-03
update): Eq. 24 is a reasonable variance proxy away from candidate zero-crossings, where
FastMspec's memory-decoupling argument stays fully intact, but should not be trusted near
candidate crossings without a literal per-crossing jackknife or a genuinely cross-spectral
closed form (Walden 2000's Wishart-based multivariate treatment is the natural place to check
next, not Haley & Anitescu's univariate theory).

**Scope of this check**: one real pair (SKRH-BAND), three `K` values. Not yet checked across the
dataset's real distance/coherence range -- flagged as an open item in the tex doc, not treated
as fully resolved by this one test.
