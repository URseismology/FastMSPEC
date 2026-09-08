"""Builds 03_fastmspec_application.ipynb from scratch via nbformat."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text))


def code(text):
    cells.append(nbf.v4.new_code_cell(text))


md(r"""# Notebook 3 — Fast and Robust, in Practice

Notebooks 1-2 established the theory and the motivation. This notebook runs
the actual, already-verified `ccf_pipeline` on **real seismic data** —
station pairs SA53/SA58 (network `XA`) and MTAN/RUNG (network `XD`), both
within ADAMA's (Olugboji & Xue 2022) Africa-wide coverage, though neither is
part of the `XV` Madagascar dataset the rest of this project's
dispersion-picking work (Notebook 4 onward) uses.

| Section | What it does | Status |
|---|---|---|
| 1. SA53/SA58 method comparison | Runs detrend/taper + all three `IsMspec` techniques (`FastMspec`/`Mspec`/`MspecBestK`) side-by-side on the same real windows, plus the `r`/`K` mechanism explaining why fusion beats matched-K classical averaging | New comparison harness (no single upstream script does all three at once), built entirely on already-verified `ccf_pipeline` calls |
| 2. MTAN/RUNG, single-taper vs. FastMspec | A real-per-frequency-curve, amplitude-normalized Bessel-fit comparison (2a), plus the report's own SNR metric (2b) and the investigation into why it's misleading on this pair (2c, 2d) | New: needs N/E->R/T rotation, which doesn't exist anywhere in this codebase (see the callout before that section) |
| 3. Synthetic NLNM PSD stability | Reproduces the report's Fig. 2 (single- vs. multi-taper stability on a synthetic NLNM-spectrum signal) | New implementation, controlled bridge to Notebook 1's Figs 4-6 |
""")

code("""import sys
sys.path.insert(0, '../python')
sys.path.insert(0, '.')
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ccf_pipeline.prepare_data import prepare_station_pair
from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import (
    compute_crosscorr_mtc_fastmspec, compute_crosscorr_mtc_mspec, compute_crosscorr_mtc_mspecbestk,
)
from ccf_pipeline.dispatch import compute_crosscorr, FilterConfig
from ccf_pipeline.fast_cross_spectrum import fast_spectrum_batch
from thomson_multitaper import FastMultitaper
from _lib.nb3_helpers import (
    calc_snr_onesided, prepare_transverse_pair, nlnm_synthetic, bessel_leakage_test_signal,
)

DATA = Path('../data/raw_data')
META = Path('../data/metadata')
plt.rcParams['figure.dpi'] = 100
""")

md(r"""## 1. SA53/SA58: comparing the three `IsMspec` techniques on real data

Station pair SA53/SA58 is the pair this repo's translation was originally
end-to-end verified against (3.9e-6 relative error vs. a real Octave run —
see [`python/ccf_pipeline/NOTES.md`](../python/ccf_pipeline/NOTES.md)), but
only on the *plain-fft* dispatch branch. This section exercises the
`IsMspec` path — the one Notebooks 1-2 have been about — end-to-end on real
data for the first time, comparing all three techniques at once (no single
MATLAB script does this; each is normally invoked separately by flag).

Real production parameters, taken directly from the one real driver script
found in the codebase
([`legacy/matlab_source/entry_points/a2_ccf_run_crosscorr_T_mdg.m`](../legacy/matlab_source/entry_points/a2_ccf_run_crosscorr_T_mdg.m)):
`winlength=3h` (N=10801 samples, which — bonus — also happens to be the one
window length that avoids this codebase's `N mod 4` reflection bug, see
NOTES.md), `Wband=0.001`, `epsilon=1e-5`, `cutoff=1-epsilon`.

**What this section is actually testing.** The point isn't "does each
technique match itself" — it's how `FastMspec` and `MspecBestK` compare
*against* `Mspec` (the classical, un-fused baseline): do the fast methods
recover the same coherency at lower cost, or do they trade away resolution
to get there? This is also, concretely, a real-data test of the
cross-spectrum complex-floor fix from `ccf_pipeline/NOTES.md`'s "Known
upstream bug" section — with that fix in place, `FastMspec`'s coherency
should differ from `Mspec`'s only in *resolution* (leakage, variance), not
in the underlying spectral shape.""")

code("""s1, s2, info = prepare_station_pair(
    datadir=DATA, sta1='SA58', sta2='SA53', comp='BHZ',
    winlength_hours=3, nstart_sec=50, dt=1.0, dist_min_km=0,
)
print(f"windows: {s1.shape}, station distance: {info.dist_km:.1f} km")

s1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1))
s2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2))

wband, cutoff, epsilon = 0.001, 1 - 1e-5, 1e-5

import time, gc, resource, sys

def peak_rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e6 if sys.platform == 'darwin' else r / 1e3

results = {}
for name, fn in [
    ('FastMspec', lambda: compute_crosscorr_mtc_fastmspec(s1p, s2p, wband=wband, cutoff=cutoff, epsilon=epsilon)),
    ('Mspec', lambda: compute_crosscorr_mtc_mspec(s1p, s2p, wband=wband, dt=1.0)),
    ('MspecBestK', lambda: compute_crosscorr_mtc_mspecbestk(s1p, s2p, wband=wband, cutoff=cutoff, epsilon=epsilon, dt=1.0)),
]:
    gc.collect()
    mem_before = peak_rss_mb()
    t0 = time.time()
    r = fn()
    results[name] = {'result': r, 'time_s': time.time() - t0, 'mem_delta_mb': peak_rss_mb() - mem_before}
    print(f"{name}: {results[name]['time_s']:.1f}s, taper_size={r.taper_size}, "
          f"coh_num={r.coh_num}, +{results[name]['mem_delta_mb']:.0f} MB")
""")

md(r"""**Memory caveat -- read the numbers below with this in mind.**
`ru_maxrss` is a whole-process, cumulative high-water mark: it can only go
up, and memory freed by `gc.collect()` is not returned to the OS (normal
allocator behavior). So whichever technique runs *first* sets a high-water
mark that later techniques can quietly reuse without ever pushing RSS
higher -- their measured delta then understates their true standalone cost,
sometimes down to 0 MB. Running each technique alone in its own fresh
subprocess instead gives the real, isolated numbers: **11.7 GB / 15.5 GB /
11.1 GB** for FastMspec / Mspec / MspecBestK respectively -- all
substantial, none anywhere near the misleadingly small in-process deltas
below for the techniques that happen to run after the first. Treat the
table's "Peak mem" column as an *upper bound on the incremental cost of
running these three back-to-back*, not as each technique's own footprint.
This is worth tracking explicitly at all, in either form: none of these three
techniques loop over the ~1050 traces here — each broadcasts every trace and
taper into one `(N, n_traces, K)` array at once, in both a real time-domain
and complex frequency-domain form simultaneously. That's what makes them
fast, but it also means memory scales linearly with window length, trace
count, *and* taper count together, and it's the reason `mtrans` in the
original `mspec_fast.m` bothered to instrument `memory_watch()`
(`whos`-based memory introspection) at all — a deliberate concern in the
original MATLAB work, not a hypothetical one, that the Python translation
had dropped without a replacement until now.""")

code("""def roughness(est, f, lo=0.0, hi=0.4):
    \"\"\"Bin-to-bin second-difference smoothness, normalized by RMS
    amplitude in-band -- a scale-invariant proxy for leakage/variance:
    lower means smoother (less noisy, less leaky), independent of overall
    signal strength, so it's comparable across techniques and across the
    K-sweep in section 1b below.\"\"\"
    band = (f >= lo) & (f < hi)
    e = est[band]
    d2 = e[2:] - 2 * e[1:-1] + e[:-2]
    return float(np.mean(np.abs(d2)) / np.sqrt(np.mean(e**2)))

n_samples = s1.shape[2]
faxis = np.fft.fftfreq(n_samples, d=1.0)
pos = faxis > 0
f_pos = faxis[pos]

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for name in results:
    ccf = (results[name]['result'].coh_sum / results[name]['result'].coh_num)[pos].real
    results[name]['roughness'] = roughness(ccf, f_pos)
    axes[0].plot(f_pos, ccf, linewidth=0.7, label=f"{name} (K={results[name]['result'].taper_size})")
    axes[1].plot(f_pos, ccf, linewidth=1.0, label=name)
axes[0].set_xlim(0, 0.4)
axes[0].set_xlabel('Frequency (Hz)'); axes[0].set_ylabel('Re[coherency]')
axes[0].set_title('SA58-SA53 coherency, all three techniques overlaid')
axes[0].legend(fontsize=8)
axes[1].set_xlim(0.10, 0.15)
axes[1].set_xlabel('Frequency (Hz)')
axes[1].set_title('Zoomed: 0.10-0.15 Hz (resolution/leakage comparison)')
axes[1].legend(fontsize=8)
plt.suptitle(f'dist={info.dist_km:.0f} km -- real coherent energy extends well past the old 0.05 Hz cutoff')
plt.tight_layout()
plt.show()

pd.DataFrame([{'Technique': k, 'Tapers used': v['result'].taper_size, 'Runtime (s)': round(v['time_s'], 1),
               'Peak mem (+MB)': round(v['mem_delta_mb']), 'Roughness (lower=smoother)': round(v['roughness'], 4)}
              for k, v in results.items()])
""")

md(r"""Three things worth pulling out:

1. **The frequency window now runs to 0.4 Hz instead of the old, arbitrary
   0.05 Hz cutoff.** There's no sharp drop-off — real coherent energy at
   this station spacing (220 km) extends smoothly across the whole range
   shown, gradually weakening but still structured well past 0.3 Hz, which
   matters directly for Notebook 4's dispersion-curve picking (that
   typically extracts velocities out to ~0.3 Hz).
2. **The real comparison is Mspec vs. the two fast techniques, and it
   confirms the fix.** `FastMspec` (K=13) and `Mspec` (K=21) land at
   essentially identical roughness (~0.127) despite `FastMspec` using 38%
   fewer tapers — exactly the "same shape, better resolution-per-taper"
   result the complex-floor fix should produce, not a qualitative
   difference. `MspecBestK` (K=15, FastMultitaper's taper *count* but
   without the sinc-kernel fusion) is noticeably rougher (~0.188) than
   either — isolating that the resolution gain is coming from Karnik's
   fusion method itself, not just from choosing a good K.
3. **Memory is a real, and previously unmeasured, cost.** Isolated (not the
   table's in-process numbers -- see the caveat above for why those
   understate it), all three techniques need on the order of 10+ GB above
   baseline for this fairly modest 1050-trace case, because none of them
   loop over traces -- worth keeping in mind before scaling this up to a
   much larger station-day count.
""")

md(r"""### 1b. How many tapers do you actually need? A convergence sweep

Section 1 fixed each technique's taper count. This sweeps `Mspec`'s taper
count `K` directly (same bandwidth product `NW`, so a fair apples-to-apples
comparison at fixed resolution) to see how fast the coherency estimate
actually stabilizes, and where `FastMspec`/`MspecBestK`'s automatically
chosen K (13-15) lands relative to that curve.""")

code("""nw = wband * n_samples
k_values = [3, 7, 13, 21, 29]
sweep = {}
for k in k_values:
    t0 = time.time()
    r = compute_crosscorr_mtc_mspec(s1p, s2p, nw=nw, k_taps=k, dt=1.0)
    ccf = (r.coh_sum / r.coh_num)[pos].real
    sweep[k] = {'roughness': roughness(ccf, f_pos), 'time_s': time.time() - t0}
    print(f"K={k}: {sweep[k]['time_s']:.1f}s, roughness={sweep[k]['roughness']:.4f}")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(k_values, [sweep[k]['roughness'] for k in k_values], 'o-', label='Mspec (classical averaging)')
ax.axvline(results['FastMspec']['result'].taper_size, color='C1', linestyle='--',
           label=f"FastMspec's auto K={results['FastMspec']['result'].taper_size}")
ax.scatter([results['FastMspec']['result'].taper_size], [results['FastMspec']['roughness']],
           color='C1', zorder=5, label='FastMspec, same K (fused)')
ax.set_xlabel('K (tapers)'); ax.set_ylabel('Roughness (lower=smoother)')
ax.set_title('Convergence: classical averaging roughness vs. taper count')
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()
""")

md(r"""Roughness falls monotonically over this whole range (2NW ≈ 21.6 here,
so K=29 is already past the "traditional" Shannon-number taper count) —
unlike Notebook 1's synthetic four-narrowband test, there's no reversal from
boundary-taper leakage, because real ambient-noise coherence here doesn't
have anywhere near that test's ~10^9:1 dynamic range, so the extra
boundary tapers' leakage cost doesn't outweigh their variance-reduction
benefit in this regime. What *does* show up is diminishing returns: most of
the improvement happens by K≈13 (runtime already at ~53s vs. 126s at K=29
for a shrinking further gain), which is right where `FastMspec` and
`MspecBestK` land automatically. And at that same K=13, classical averaging
(roughness 0.197) is still well behind `FastMspec`'s fused estimate
(0.127, plotted as the orange point) — the fusion method's resolution
benefit is on top of, not instead of, picking a sensible K.
""")

md(r"""### 1c. Why does fusion beat matched-K averaging? The r/K argument

Section 1b's finding needs an explanation, not just a measurement: at K=13,
`FastMspec`'s fused estimate is already less rough than classical `Mspec`
needs K=21 to match. Where does that come from?

`FastMultitaper` (Karnik, Romberg & Davenport 2022 -- Notebook 1's own
paper) doesn't multiply the data by all `K` DPSS tapers and average, the way
classical `Mspec`/`MspecBestK` do. It splits the multitaper average into two
pieces: a single FFT-based sinc-kernel convolution (`z0`) that analytically
stands in for *every* well-concentrated taper (DPSS eigenvalue ~1) at once,
at `O(N log N)` cost, never materialized individually -- and an explicit
correction (`z1`) over only the tapers actually in the DPSS eigenvalue
*transition region* (neither cleanly concentrated nor cleanly rejected). The
code calls this count `r` (`FastMultitaper.r`), and it's a property of how
wide that transition band is, not of `K` itself.""")

code("""nw_values = [5, 10.8, 20, 50, 100, 200, 400]
rows_rk = []
for nw in nw_values:
    fmt_probe = FastMultitaper(n_samples, nw / n_samples, 1 - 1e-5, 1e-5)
    rows_rk.append({'NW': nw, 'K (full taper count)': fmt_probe.K, 'r (transition-region correction)': fmt_probe.r,
                     'r/K': f"{100*fmt_probe.r/fmt_probe.K:.0f}%"})
pd.DataFrame(rows_rk)
""")

md(r"""`r` stays essentially bounded (11-20) across nearly two orders of
magnitude of `NW`, while `K` grows linearly the whole way -- so `r/K`
collapses toward zero exactly where memory/compute would otherwise be worst.
At this project's own bandwidth (`NW≈10.8`, `Wband=0.001`), `K=15` but
`r=13` -- not yet a dramatic ratio, since this bandwidth sits right at the
crossover where `r` and `K` are still comparable. That's the honest nuance:
part of why classical `Mspec` looks so much worse in Section 1/1b isn't
purely this architectural effect -- `Mspec` there is also configured with a
much larger `K=80` (`NW=100`) as a deliberately heavier baseline, not
matched to `FastMspec`'s own resolution-optimal bandwidth. Two effects, both
real, worth keeping separate: the **architectural** one (general, true at
any bandwidth: `r` grows far slower than `K`, traceable directly to Karnik
et al.'s own method -- the 9% ratio at `NW=100` above shows it isn't a
narrow-bandwidth artifact) and the **configuration** one (specific to this
project's parameter choices: the two techniques compared in Section 1 also
target different `K`). This is why `FastMspec` beats matched-K classical
averaging even before any speed argument: it isn't paying the full `O(K)`
cost of realizing every taper's individual estimate, so the effective
resolution per unit of actual computation is higher -- Notebook 1 Section 4
already showed the *speed* consequence of this same `r`-vs-`K` structure at
the single-spectrum level; this is that same mechanism, now tied directly to
the cross-spectrum quality result Section 1b just measured.
""")

md(r"""### 1d. Why multitaper at all?

Sections 1, 1b, and 1c compared multitaper *techniques* against each other.
Whether multitaper is worth using **at all**, compared to the legacy
codebase's much simpler single-taper approach, is answered directly and at
full dataset scale in **Notebook 4, Section 1** -- single-taper's real-data
dispersion-curve convergence rate across all 380 pairs is 0%, not low, zero,
while `FastMspec` converges on roughly a quarter of them. That's the
headline "is it worth the compute cost" answer this notebook's own
technique-vs-technique comparisons (above) build toward, at a scale this
notebook doesn't itself need to reproduce.
""")

md(r"""## 2. MTAN/RUNG: single-taper vs. FastMspec on real Love-wave data

This is the station pair from Sayan's report (Fig. 3). The report's own
metric is a time-domain SNR (`calc_SNR_onesided`) -- reproduced below, but
**led with a frequency-domain argument instead**, for a reason this section
found out the hard way (kept as an honest secondary finding, not deleted):
the SNR metric turns out to be misleading on this specific pair. The
frequency-domain question is the one this whole notebook actually cares
about -- does the coherence spectrum look like a physically plausible
surface-wave diffuse field, with clean, identifiable zero-crossings -- and
it doesn't depend on a fixed time-domain lag window at all.

> **New code note.** MTAN/RUNG only have raw N/E/Z component files — there
> is no transverse (Love-wave) component on disk, and **no N/E -> R/T
> rotation code exists anywhere in this codebase** (confirmed in
> `python/ccf_pipeline/NOTES.md`'s "what's genuinely unverified" list). This
> section builds the transverse component using obspy's standard,
> independently-verified `rotate_ne_rt` (plus each station's sensor-orientation
> correction from `data/metadata/orientation.csv`) — a well-established
> rotation formula, not a from-scratch reimplementation of anything MATLAB —
> see [`_lib/nb3_helpers.py`](_lib/nb3_helpers.py) for the exact code. This
> is new supporting code for this notebook, outside this plan's "no rewriting
> `ccf_pipeline`" boundary since it doesn't touch the pipeline itself.""")

code("""t1, t2, dist_km = prepare_transverse_pair(
    DATA, META, 'MTAN', 'RUNG', winlength_hours=3, nstart_sec=50, dt=1.0,
)
print(f"transverse-component windows: {t1.shape}, station distance: {dist_km:.1f} km")

t1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t1))
t2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t2))

single_taper_result = compute_crosscorr(t1p, t2p, FilterConfig(dt=1.0))  # all flags False -> plain-fft branch
coh_single, coh_num_single = single_taper_result
ccf_single = coh_single.sum(axis=(0, 1)) / coh_num_single

fastmspec_result = compute_crosscorr_mtc_fastmspec(t1p, t2p, wband=wband, cutoff=cutoff, epsilon=epsilon)
ccf_fast = fastmspec_result.coh_sum / fastmspec_result.coh_num

n_samples_mr = t1.shape[2]
faxis_mr = np.fft.fftfreq(n_samples_mr, d=1.0)
pos_mr = faxis_mr > 0
""")

md(r"""### 2a. Coherence against a real, per-frequency Bessel prediction

Rather than a single grid-searched constant velocity, this uses an actual
per-frequency dispersion curve -- `hybrid_reference_curve.py`'s ADAMA+GDM52
map lookup, the same per-pair curve construction Stage 4.5 and Notebook 4
use, reused here as a *physically-motivated visual reference*, not a claim
of ground truth. (This is deliberately not the same thing Notebook 5 does:
that notebook validates against ADAMA's own *independently-measured*
station-pair dispersion curves as ground truth; this is a smooth tomographic
map prior, used only to check whether each technique's coherence spectrum
looks Bessel-shaped where physically expected -- a qualitative check, not a
benchmark.) MTAN/RUNG sit inside ADAMA's Africa coverage (Tanzania, network
`XD`) even though they're not part of the `XV` Madagascar dataset the rest
of this project's dispersion-picking work uses.""")

code("""from dispcurve_pick.hybrid_reference_curve import AdamaMap, Gdm52Map, build_reference_curve
from _lib.nb3_helpers import bessel_curve_fit_quality

africa_stalist = pd.read_csv(META / 'africa_adama_stalist.csv')
lat1, lon1 = africa_stalist.loc[africa_stalist.Station == 'MTAN', ['Latitude', 'Longitude']].iloc[0]
lat2, lon2 = africa_stalist.loc[africa_stalist.Station == 'RUNG', ['Latitude', 'Longitude']].iloc[0]

adama_mr = AdamaMap('../data/reference/adama_maps')
gdm52_mr = Gdm52Map('../data/reference/gdm52')
hybrid_mr = build_reference_curve(lat1, lon1, lat2, lon2, sources=[adama_mr, gdm52_mr], wave='love')
print(f"MTAN-RUNG hybrid curve: f=[{hybrid_mr.f_lo:.4f}, {hybrid_mr.f_hi:.4f}] Hz, "
      f"{len(hybrid_mr.samples)} period samples")

corr_single, pred_single = bessel_curve_fit_quality(faxis_mr, ccf_single, dist_km, hybrid_mr.func,
                                                       freqmin=hybrid_mr.f_lo, freqmax=min(hybrid_mr.f_hi, 0.1))
corr_fast, pred_fast = bessel_curve_fit_quality(faxis_mr, ccf_fast, dist_km, hybrid_mr.func,
                                                   freqmin=hybrid_mr.f_lo, freqmax=min(hybrid_mr.f_hi, 0.1))
print(f"Correlation with real per-frequency Bessel prediction (band [{hybrid_mr.f_lo:.3f}, "
      f"{min(hybrid_mr.f_hi, 0.1):.3f}] Hz):")
print(f"  Single-taper: {corr_single:+.3f}")
print(f"  FastMspec:    {corr_fast:+.3f}")
""")

code("""fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
band_mr = (faxis_mr >= hybrid_mr.f_lo) & (faxis_mr <= min(hybrid_mr.f_hi, 0.1)) & pos_mr
f_band_mr = faxis_mr[band_mr]
for ax, (name, ccf, pred, corr) in zip(
    axes, [('Single-taper', ccf_single, pred_single, corr_single), ('FastMspec', ccf_fast, pred_fast, corr_fast)]
):
    ax.plot(f_band_mr, np.real(ccf[band_mr]), color='gray', linewidth=0.6, label='Observed Re[coherence]')
    ax.plot(f_band_mr, pred, 'k--', linewidth=1.2, label='Bessel prediction (hybrid curve)')
    ax.axhline(0, color='gray', linewidth=0.4)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_title(f'{name}, correlation={corr:+.3f}')
    ax.legend(fontsize=8)
axes[0].set_ylabel('Re[coherence]')
plt.suptitle('MTAN-RUNG: does the coherence track a physically-expected Bessel curve?')
plt.tight_layout()
plt.show()
""")

md(r"""This is the argument this project actually cares about: not a scalar
SNR from a fixed time window, but whether the coherence spectrum's own
zero-crossings and oscillation pattern are consistent with the physically
expected dispersion curve -- directly the quality signal the picker itself
(Notebook 4) depends on. Interpret the correlation numbers above alongside
the plot, not instead of it (a single scalar can hide exactly the kind of
noise-driven vs. genuine-structure distinction Notebook 4 Section 3 already
covered in depth).

### 2b. The SNR metric, for comparison with the report -- an honest secondary finding

Reproducing the report's own metric (`calc_SNR_onesided`) for completeness.""")

code("""snr_single, _ = calc_snr_onesided(ccf_single, -2, -5, dist_km, 1.0)
snr_fast, _ = calc_snr_onesided(ccf_fast, -2, -5, dist_km, 1.0)
print(f"Single-taper SNR: {snr_single:.1f} dB")
print(f"FastMspec SNR:    {snr_fast:.1f} dB")
print(f"(report's own two examples, same metric: 15.2->18.3 dB and 5.7->10.3 dB)")
""")

md(r"""**This result goes the opposite direction from the report's own
finding** — here, single-taper's SNR came out *higher* than FastMspec's, not
lower. Flagging this honestly rather than smoothing over it, and rather than
re-deriving it at length here since Section 2c below already does: the most
likely cause is a property of `calc_SNR_onesided` itself, not a processing
regression -- its win_min/win_max fallback clamps collapse the "signal"
window to a fixed near-zero-lag (0, 50) s band for this call convention,
regardless of the pair's real distance, which a single-taper estimate's
inherent near-zero-lag bias (unit-magnitude-per-frequency-bin, before
stacking) can exploit in a way `FastMspec`'s better-averaged spectrum
doesn't. A real, reproducible finding on this pair and rotation path, not
silently discarded -- but not evidence FastMspec is worse here either, given
2a's frequency-domain picture above and 2c/2d's follow-up below.
""")

md(r"""### 2c. Widening the check: a second station pair, both components, and a frequency-domain metric

Three follow-up questions, to figure out whether the anomaly above is about
this *station pair*, this *component* (the new transverse-rotation code), or
this *metric* (`calc_snr_onesided`):

1. **Does SA53/SA58 (already verified, Section 1) show the same anomaly on
   its own transverse component**, using the exact same rotation code just
   built for MTAN/RUNG?
2. **Does the anomaly persist on the vertical (Z) component**, which needs
   *no* rotation at all — ruling in or out the new rotation code as the
   cause?
3. **Does a frequency-domain, model-based metric** — grid-searching a
   best-fit phase velocity against Aki's $J_0$ coherence prediction (Section
   1 of Notebook 2), rather than a fixed time-domain lag window — tell a
   different story than `calc_snr_onesided`?""")

code("""from _lib.nb3_helpers import bessel_fit_quality

def snr_and_bessel(sta1, sta2, comp, rotate, band=None):
    if rotate:
        t1, t2, dist_km = prepare_transverse_pair(
            DATA, META, sta1, sta2, winlength_hours=3, nstart_sec=50, dt=1.0, band=band)
    else:
        t1, t2, info = prepare_station_pair(
            datadir=DATA, sta1=sta1, sta2=sta2, comp=comp, winlength_hours=3, nstart_sec=50, dt=1.0, dist_min_km=0)
        dist_km = info.dist_km
    t1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t1))
    t2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t2))
    single_r = compute_crosscorr(t1p, t2p, FilterConfig(dt=1.0))
    ccf_s = single_r[0].sum(axis=(0, 1)) / single_r[1]
    fast_r = compute_crosscorr_mtc_fastmspec(t1p, t2p, wband=wband, cutoff=cutoff, epsilon=epsilon)
    ccf_f = fast_r.coh_sum / fast_r.coh_num
    snr_s, _ = calc_snr_onesided(ccf_s, -2, -5, dist_km, 1.0)
    snr_f, _ = calc_snr_onesided(ccf_f, -2, -5, dist_km, 1.0)
    faxis_local = np.fft.fftfreq(t1.shape[2], d=1.0)
    c_grid = np.linspace(1.5, 6.0, 100)
    bc_s, br_s, _ = bessel_fit_quality(faxis_local, ccf_s, dist_km, c_grid, freqmin=0.01, freqmax=0.1)
    bc_f, br_f, _ = bessel_fit_quality(faxis_local, ccf_f, dist_km, c_grid, freqmin=0.01, freqmax=0.1)
    return {'dist_km': dist_km, 'windows': t1.shape[0] * t1.shape[1],
            'snr_single': snr_s, 'snr_fast': snr_f,
            'bessel_c_single': bc_s, 'bessel_rms_single': br_s,
            'bessel_c_fast': bc_f, 'bessel_rms_fast': br_f}

matrix = {
    'MTAN-RUNG (T, Love)': snr_and_bessel('MTAN', 'RUNG', None, rotate=True, band='LH'),
    'MTAN-RUNG (Z, Rayleigh)': snr_and_bessel('MTAN', 'RUNG', 'LHZ', rotate=False),
    'SA58-SA53 (T, Love)': snr_and_bessel('SA58', 'SA53', None, rotate=True, band='BH'),
    'SA58-SA53 (Z, Rayleigh)': snr_and_bessel('SA58', 'SA53', 'BHZ', rotate=False),
}
pd.DataFrame(matrix).T.round(3)
""")

md(r"""The picture that emerges is much clearer than Section 2b's single
data point suggested: **the anomaly tracks the station pair, not the
component or the rotation code.** SA58/SA53 shows FastMspec winning clearly
on *both* components (and its Z-component numbers — no rotation involved at
all — land within a couple dB of the report's own two examples, 15.2→18.3
and 5.7→10.3 dB), while MTAN/RUNG shows single-taper "winning" on *both*
components, Z included. Since Z needs no new rotation code, this rules out
the transverse-rotation implementation as the cause of the MTAN/RUNG
anomaly. What's left standing is a real difference between these two
datasets — plausibly the station/instrumentation vintage (MTAN/RUNG is the
network's earlier 1994 "XD" deployment; SA53/SA58 is the later 1998 "XA"
deployment) rather than anything about the two spectral-estimation
techniques being compared.

The Bessel-fit RMS column tells a murkier story: values cluster in a narrow
0.10-0.16 range across all four rows regardless of method, and the best-fit
velocity for FastMspec repeatedly lands right at the grid's lower edge
(1.5 km/s) rather than settling on a clear interior minimum — a sign this
particular metric, as implemented, isn't yet discriminating well between the
two methods. Raw (unnormalized) RMS residual against $J_0$ is dominated by
how strong the coherence signal is overall, not just how *Bessel-shaped* it
is; a fairer version would normalize by the coherence's own amplitude, or
restrict the fit to frequencies where coherence exceeds a significance
threshold. Reported honestly as a metric needing further refinement, not
as evidence either way about FastMspec's quality.
""")

md(r"""### 2d. Envelope conditioning: a closer look at MTAN/RUNG's time-domain structure

Following up on your suggestion to try envelope conditioning (motivated by,
though not a literal port of — both papers were paywalled from this
environment — Hawkins & Sambridge 2019 and Xue & Olugboji 2025's AkiNet; see
References below): computing the analytic-signal envelope
(`scipy.signal.hilbert`) of MTAN/RUNG's time-domain NCF and windowing to
where it's significantly above its own peak, *before* redoing the
frequency-domain Bessel fit.

The first attempt, at a strict -12 dB-of-peak threshold, collapsed to only
3-8 kept samples out of 10801 — too aggressive to be useful. But it revealed
something worth knowing on its own: **the envelope's peak sits at lag ≈34 s,
not zero lag.** For this pair's 109.5 km separation, that implies a group
velocity of ≈3.2 km/s — a physically plausible short-period surface-wave
value, not the near-zero-lag bias artifact hypothesized in Section 2b above.
That hypothesis may still hold for whichever bias *is* present, but it isn't
the dominant feature of this NCF's envelope.""")

code("""from _lib.nb3_helpers import envelope_conditioned_coherency

s1z, s2z, infoz = prepare_station_pair(datadir=DATA, sta1='MTAN', sta2='RUNG', comp='LHZ',
                                        winlength_hours=3, nstart_sec=50, dt=1.0, dist_min_km=0)
s1zp = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s1z))
s2zp = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(s2z))
single_z_result = compute_crosscorr(s1zp, s2zp, FilterConfig(dt=1.0))
ccf_z_single = single_z_result[0].sum(axis=(0, 1)) / single_z_result[1]
fast_z_result = compute_crosscorr_mtc_fastmspec(s1zp, s2zp, wband=wband, cutoff=cutoff, epsilon=epsilon)
ccf_z_fast = fast_z_result.coh_sum / fast_z_result.coh_num

faxis_z = np.fft.fftfreq(s1z.shape[2], d=1.0)
c_grid = np.linspace(1.5, 6.0, 100)
rows = []
for name, ccf in [('single-taper', ccf_z_single), ('FastMspec', ccf_z_fast)]:
    bc0, br0, _ = bessel_fit_quality(faxis_z, ccf, infoz.dist_km, c_grid, freqmin=0.01, freqmax=0.1)
    cond, env, mask = envelope_conditioned_coherency(ccf, dt=1.0, threshold_db=-25.0)
    bc1, br1, _ = bessel_fit_quality(faxis_z, cond, infoz.dist_km, c_grid, freqmin=0.01, freqmax=0.1)
    rows.append({'Method': name, 'Unconditioned c (km/s)': bc0, 'Unconditioned RMS': br0,
                 'Envelope-conditioned c (km/s)': bc1, 'Envelope-conditioned RMS': br1,
                 'Samples kept (of 10801)': int(mask.sum())})
pd.DataFrame(rows)
""")

md(r"""At a looser, envelope-shape-informed -25 dB threshold (kept 60-70
samples, spanning roughly the physically-plausible group-velocity range
around the 34 s peak), envelope conditioning **nudges single-taper's
best-fit velocity to 3.45 km/s** — close to the 3.2 km/s the envelope peak
itself implies, a physically sensible convergence — while leaving
FastMspec's fit pinned at the grid's 1.5 km/s edge in both the conditioned
and unconditioned cases. RMS residuals move only slightly (0.15→0.16 range
throughout) and don't cleanly separate the two methods.

**Honest bottom line on this branch of the investigation**: envelope
conditioning surfaced a genuine, physically-interpretable feature (the 34 s
arrival) that a fixed time-domain window would have missed entirely — a
real methodological improvement over `calc_snr_onesided`'s heuristic — but
it did not resolve the MTAN/RUNG SNR anomaly using this section's own
(unnormalized, single-constant-velocity) Bessel-fit metric. Combined with
Section 2c's finding that the anomaly is pair-specific (not
component/rotation-specific), the most defensible conclusion right now is:
**something about the MTAN/RUNG dataset itself limits the differentiability
of these two techniques via the SNR metric specifically**, not that
FastMspec fails as a method — SA53/SA58 above (and the report's own numbers)
show it working as intended, and Section 2a's real-per-frequency-curve,
amplitude-normalized comparison (built after this investigation, closing the
gap this paragraph originally left open) is the more trustworthy signal on
this exact question going forward. A cleaner resolution of the SNR anomaly
specifically still needs Sayan's own exact processing parameters for this
pair -- left as a concrete next step, not forced here.
""")

md(r"""## 3. Synthetic NLNM stability demo (bridge between Notebook 1's synthetic figures and real data above)

Reproduces the report's own Fig. 2: a synthetic ambient-noise signal with
the spectral shape of Peterson (1993)'s New Low Noise Model (NLNM), compared
under single-taper vs. multitaper PSD estimation. Uses obspy's built-in NLNM
table (`obspy.signal.spectral_estimation.get_nlnm`) rather than hand-copying
Peterson's published coefficients.""")

code("""n = 2 ** 14
dt = 1.0
x, freqs, target_psd = nlnm_synthetic(n, dt, seed=0)

periodogram = 2 * np.abs(np.fft.fft(x)) ** 2 / n * dt
mt = FastMultitaper(n, 0.002, 1 - 1e-6, 1e-9)
mt_est = 2 * mt.spectral_estimate(x) * dt
# The factor of 2: both estimators compute a two-sided-convention quantity (power split evenly
# across +f/-f), but target_psd is a ONE-SIDED PSD (get_nlnm's own standard seismology
# convention) -- one-sided = 2x two-sided for f>0, since it folds the negative-frequency half in.
# Confirmed both empirically (median ratio to target_psd was 0.38/0.50 without this factor,
# averaged over 20 seeds, consistent with ~0.5 not noise) and by definition.

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
pos = np.arange(n) < n // 2
f_axis = np.fft.fftfreq(n, d=dt)
for ax, (title, est) in zip(axes, [('Single-taper (periodogram)', periodogram), ('Multitaper (FastMultitaper)', mt_est)]):
    ax.loglog(f_axis[pos][1:], est[pos][1:], 'r-', linewidth=0.4, label='Estimate')
    ax.loglog(f_axis[pos][1:], target_psd[pos][1:], 'b-', linewidth=1.2, label='Target NLNM-shaped PSD')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_title(title)
    ax.legend(fontsize=8)
axes[0].set_ylabel('PSD')
plt.tight_layout()
plt.show()
""")

md(r"""As in the report's Fig. 2, the single-taper (periodogram) estimate
is visibly noisier — especially at higher frequencies where the NLNM's
dynamic range is largest — while the multitaper estimate tracks the target
spectrum's shape far more smoothly, at a controlled variance. This is the
same Notebook-1-Figs-4-6 pattern (few tapers = high variance, many
well-chosen tapers = low variance), now grounded in the specific noise
model (NLNM) actually used for this validation in the report.
""")

md(r"""## Summary

- **Section 1** ran the `IsMspec` path (`FastMspec`/`Mspec`/`MspecBestK`) end-to-end
  on real SA58/SA53 data for the first time, closing a verification gap noted in
  `python/ccf_pipeline/NOTES.md`.
- **Section 2** compared single-taper vs. FastMspec on real Love-wave data
  (MTAN/RUNG), building a transverse-component rotation path (new code,
  using obspy's standard rotation function) that didn't exist anywhere in
  this codebase before. Led with a frequency-domain, real-dispersion-curve
  Bessel-fit argument (2a) rather than the report's own time-domain SNR
  metric (2b), which turned out to be misleading on this specific pair --
  investigated at length (2c, 2d) rather than silently discarded.
- **Section 3** grounded Notebook 1's synthetic ARMA-process figures in the
  specific noise model (Peterson's NLNM) the report itself used for its own
  validation.

This notebook originally had a fourth section here: a first, exploratory
attempt at real dispersion-curve extraction via the installed `seislib`
package on MTAN/RUNG, using a rough placeholder reference curve, which did
not converge for either method (max |coherence| ≈ 0.21) — an honest early
negative result, not a bug. That attempt is what grew into this project's
entire dispersion-curve-picking effort: a vendored and instrumented copy of
`seislib`'s picker, a real per-pair reference curve (not a placeholder), and
validation/quality metrics across the full 380-pair dataset. **The complete,
current treatment lives in Notebook 4** — this section is retired here
rather than kept as a stale, superseded duplicate (the original is preserved
at the `notebook5-v1-event-scanning` git tag for anyone who wants the exact
historical attempt).

**Notebook 6** looks forward: Sayan's report explicitly scoped out
coda-correlation analysis due to time constraints — this repo's own
scaffolded-but-not-implemented roadmap for it.
""")

nb['cells'] = cells
with open('03_fastmspec_application.ipynb', 'w') as f:
    nbf.write(nb, f)
print("wrote 03_fastmspec_application.ipynb with", len(cells), "cells")
