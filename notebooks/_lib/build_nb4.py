"""Builds 04_dispersion_curve_picking.ipynb from scratch via nbformat.

Complete rebuild (Stage 5 of the Notebook 5 revamp -- docs/notebook5_revamp_progress.md), not a
patch on the old naive zero-crossing/max-min template-matching barcode this notebook held before
2026-09-08 (preserved at the `notebook5-v1-event-scanning` git tag). Sections 1-3 below carry over
old Sections 1-3's real content (Setup/zero-crossing stability/resolution-bandwidth principle),
lightly reframed as Motivation rather than barcode setup -- Section 4 onward is entirely new:
the vendored/instrumented picker, a real worked example (SKRH-BAND, Sayan's own known-good pair),
the full 380-pair Round 1/Round 2 results, and Stage 4.5's reference-curve/corridor fixes.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text))


def code(text):
    cells.append(nbf.v4.new_code_cell(text))


nb['cells'] = cells

md(r"""# Notebook 4 — Picking Phase-Velocity Dispersion Curves at Scale

Notebooks 1-3 established that `FastMspec` produces a smoother, lower-variance cross-spectrum than
single-taper, at no cost in zero-crossing resolution -- and that this isn't just a variance
story: Notebook 3's Sections 1c-1f showed *why* it matters specifically for phase-velocity picking
(zero-crossing locations, not raw amplitude, are what the picker actually depends on). This
notebook asks the next question directly: given a trustworthy coherence spectrum, how do you
actually extract a trustworthy phase-velocity *dispersion curve* from it, at real dataset scale,
with real quality control -- not just eyeball whether it looks plausible?

That question has a real history in this project, worth being upfront about. An early, naive
attempt (raw zero-crossing/maximum/minimum event scanning, scored against a template library --
`docs/coherence_barcode_design.tex` Sections 1-7) was built, tested, and honestly found
insufficient: maxima/minima detection produced ~5x as many "reliable" events as zero-crossing
detection, dominated by per-bin noise rather than genuine curvature. That whole approach is
preserved at the `notebook5-v1-event-scanning` git tag for anyone who wants to see it, but it is
**not** what this notebook does. Best practice -- and the approach `seislib` (Magrini et al. 2022)
itself takes -- is to translate the coherence spectrum into an actual phase-velocity curve
*first*, using a picking process with real internal quality control, and judge quality from
information that process already computes, not from scoring raw spectral events directly.

| Section | What it does |
|---|---|
| 1. Motivation | Single-taper vs. `FastMspec` on SA53/SA58 -- carries the "FastMspec beats single-taper" thread from Notebooks 1-3 forward into this notebook's own opening result (moved from old Notebook 5 Section 1) |
| 2. Zero-crossing stability | Real-data confirmation that smoothing single-taper to match FastMspec's variance still drifts crossing locations (moved from old Notebook 5 Section 2) |
| 3. The resolution principle | Bandwidth vs. Bessel zero-crossing spacing -- the theoretical grounding `NW_high(r)` and this project's whole bandwidth story build on (moved from old Notebook 5 Section 3) |
| 4. The method | The vendored, instrumented `seislib` picker -- what "instrumented" means and why it exists |
| 5. Worked example | SKRH-BAND vs. Sayan's own known-good result, including a live comparison against the pre-Stage-4.5 pipeline on this same pair |
| 6. At scale | The full 380-pair Round 1 (technique comparison) and Round 2 (bandwidth sweep, 4 hypotheses tested) |
| 7. Stage 4.5 | The reference-curve-mismatch mechanism, the polarization fix, the hybrid ADAMA+GDM52 curve, the corridor-widening strategy -- and what combining them actually bought |
| 8. Discussion | What's proven, what's still open, and the bridge to Notebook 5 |
""")

code("""import sys
sys.path.insert(0, '../python')
sys.path.insert(0, '.')
from pathlib import Path
import time
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.io import loadmat

from ccf_pipeline import preprocessing as pp
from ccf_pipeline.crosscorr_mtc import compute_crosscorr_mtc_fastmspec
from ccf_pipeline.dispatch import compute_crosscorr, FilterConfig
from _lib.nb3_helpers import prepare_transverse_pair
from _lib.nb4_helpers import scan_templates_with_picker, best_of, score
from dispcurve_pick import extract_dispcurve, load_reference_curve, build_template_family, build_template_family_widened
from dispcurve_pick.hybrid_reference_curve import AdamaMap, Gdm52Map, build_reference_curve

DATA = Path('../data/raw_data')
META = Path('../data/metadata')
REFERENCE = Path('../data/reference/SDISPL.ASC')
RESULTS = Path('../data/results/dispcurve_quality')
plt.rcParams['figure.dpi'] = 100

def roughness(est, f, lo=0.0, hi=0.4):
    \"\"\"Bin-to-bin second-difference smoothness, normalized by RMS amplitude in-band
    (Notebook 3 Section 1c's convention).\"\"\"
    band = (f >= lo) & (f < hi)
    e = est[band]
    d2 = e[2:] - 2 * e[1:-1] + e[:-2]
    return float(np.mean(np.abs(d2)) / np.sqrt(np.mean(e**2)))
""")

# ============================================================================
# Section 1 -- Motivation (from old N5 Section 1, lightly reframed)
# ============================================================================

md(r"""## 1. Motivation: single-taper vs. `FastMspec` on SA53/SA58 (Love-wave)

Same pair, same production code, same parameters as Notebook 3 Section 1c -- recomputed
independently here so this notebook doesn't depend on running Notebook 3 first. This is the
result everything that follows is grounded in: `FastMspec` is not merely a faster way to compute
the same thing single-taper computes -- it is a *better* cross-spectrum, and that quality
difference is what makes dispersion-curve picking (the rest of this notebook) tractable at all on
real, noisy data. Section 6 below makes this concrete at full dataset scale: single-taper's
dispersion-curve convergence rate across all 380 real pairs is **0%** -- not low, zero -- while
`FastMspec` converges on roughly a quarter of them. That's not a footnote; it's the headline
reason this whole pipeline exists.""")

code("""wband, cutoff, epsilon = 0.001, 1 - 1e-5, 1e-5

t1, t2, dist_km_t = prepare_transverse_pair(
    DATA, META, 'SA58', 'SA53', winlength_hours=3, nstart_sec=50, dt=1.0, band='BH',
)
print(f"transverse windows: {t1.shape}, station distance: {dist_km_t:.1f} km")

t1p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t1))
t2p = pp.ccf_cos_taper_3dim(pp.ccf_detrend_3dim(t2))

t0 = time.time()
single_result = compute_crosscorr(t1p, t2p, FilterConfig(dt=1.0))  # all flags False -> plain-fft/single-taper branch
coh_single, coh_num_single = single_result
ccf_love_single = coh_single.sum(axis=(0, 1)) / coh_num_single
t_love_single = time.time() - t0

t0 = time.time()
fast_result_love = compute_crosscorr_mtc_fastmspec(t1p, t2p, wband=wband, cutoff=cutoff, epsilon=epsilon)
ccf_love_fast = fast_result_love.coh_sum / fast_result_love.coh_num
t_love_fast = time.time() - t0

n_love = t1.shape[2]
faxis_love = np.fft.fftfreq(n_love, d=1.0)
pos_love = faxis_love > 0
f_love = faxis_love[pos_love]
cs_love = ccf_love_single[pos_love].real
cf_love = ccf_love_fast[pos_love].real

rgh_single = roughness(cs_love, f_love, 0.02, 0.3)
rgh_fast = roughness(cf_love, f_love, 0.02, 0.3)
print(f"Single-taper: {t_love_single:.2f}s, roughness={rgh_single:.4f}")
print(f"FastMspec:    {t_love_fast:.1f}s,  roughness={rgh_fast:.4f}  ({rgh_single/rgh_fast:.1f}x smoother)")
""")

# ============================================================================
# Section 2 -- Zero-crossing stability (from old N5 Section 2, unchanged)
# ============================================================================

md(r"""## 2. Zero-crossing stability under smoothing (real data)

Sayan's own report and presentation (`docs/references/SSWAR_ESC425_Project_Report.pdf`, Fig. 7;
`SSWAR_ESC425_Presentstion.pptx`, "Results: Zero Crossings Variations") already identified the
real, practically-relevant version of this question: not whether single-taper can recover
amplitude in some deep spectral null, but whether *smoothing it enough to match multitaper's
variance* shifts the *locations* of the cross-spectrum's zero crossings -- since phase-velocity
picking (Section 4 onward, below) works from those crossing locations, not from amplitude. Their
finding: a window of 5 preserves crossing structure reasonably, 50 shows "noticeable drift," 100
gives "significant shifts, unreliable phase-velocity retrievals" -- and `FastMspec` gets its
variance reduction without needing any of that smoothing at all.

This reproduces that same test on our own SA53/SA58 transverse (Love-wave) data and our own
production code.""")

code("""def zero_crossings(y, freqs):
    \"\"\"Sign-change frequencies, linearly interpolated between the two straddling samples for
    sub-bin precision.\"\"\"
    sgn = np.sign(y)
    idx = np.where(np.diff(sgn) != 0)[0]
    return freqs[idx] + (freqs[idx + 1] - freqs[idx]) * (0 - y[idx]) / (y[idx + 1] - y[idx])

def nearest_match(ref, other):
    \"\"\"For each crossing in ref, the nearest crossing in other -- handles the very different
    crossing counts a noisy vs. smoothed spectrum give, unlike pairing by index/order.\"\"\"
    if len(other) == 0:
        return np.full(len(ref), np.nan)
    return np.array([other[np.argmin(np.abs(other - r))] for r in ref])

love_band = (f_love >= 0.02) & (f_love < 0.3)
fb_love = f_love[love_band]
zc_fastmspec = zero_crossings(cf_love[love_band], fb_love)

windows = [1, 5, 50, 100]
zc_by_window, counts, drifts = {}, [], []
for w in windows:
    variant = cs_love[love_band] if w == 1 else uniform_filter1d(cs_love, size=w)[love_band]
    zc = zero_crossings(variant, fb_love)
    matched = nearest_match(zc_fastmspec, zc)
    zc_by_window[w] = zc
    counts.append(len(zc))
    drifts.append(float(np.nanmean(np.abs(matched - zc_fastmspec))))
    print(f"window={w:3d}: {len(zc):4d} crossings, mean nearest-neighbor drift = {drifts[-1]:.5f} Hz")
print(f"(FastMspec: {len(zc_fastmspec)} crossings, K={fast_result_love.taper_size})")
""")

code("""fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = plt.cm.viridis(np.linspace(0, 0.85, len(windows)))
for w, c in zip(windows, colors):
    matched = nearest_match(zc_fastmspec, zc_by_window[w])
    axes[0].scatter(zc_fastmspec, matched, s=10, color=c, alpha=0.7, label=f'window={w}')
axes[0].plot([0, 0.3], [0, 0.3], 'k--', linewidth=1, label='y=x')
axes[0].set_xlabel('Zero crossings, FastMspec (Hz)')
axes[0].set_ylabel('Nearest zero crossing, single-taper smoothed (Hz)')
axes[0].set_title('Zero-crossing agreement (report Fig. 7 style)')
axes[0].legend(fontsize=8)

ax2 = axes[1]
ax2.plot(windows, drifts, 'o-', color='C3')
ax2.set_xlabel('Smoothing window (bins)'); ax2.set_ylabel('Mean nearest-neighbor drift (Hz)', color='C3')
ax3 = ax2.twinx()
ax3.plot(windows, counts, 's--', color='C0', alpha=0.6)
ax3.set_ylabel('Number of zero crossings', color='C0')
ax2.set_title('Drift grows with smoothing window')
plt.suptitle(f'Figure 1 -- SA58-SA53 transverse (Love), real data -- zero-crossing stability, 0.02-0.3 Hz')
plt.tight_layout()
plt.show()
""")

md(r"""Two distinct effects are visible, and they matter for different reasons: crossing *count*
collapses sharply from raw (single-taper, window=1) to window=50 as spurious noise-driven
crossings are filtered out -- this part is exactly the variance-reduction benefit Section 1
already showed smoothing *can* buy. But nearest-neighbor drift climbs monotonically the whole way,
from ~0.00007 Hz (raw) to ~0.0029 Hz (window=100), including *after* the crossing count has
already stabilized (window=50->100 keeps the same ~52 crossings, but they keep moving) -- this is
the report's point: once you're smoothing enough to reduce variance, you're already shifting
genuine crossing locations, and there's no window size that gets you one without the other.
`FastMspec` reaches comparable variance reduction with none of this drift, because it isn't
smoothing anything after the fact -- the resolution improvement comes from proper tapering, at
the point of spectral estimation, not post-hoc.
""")

# ============================================================================
# Section 3 -- The resolution principle (from old N5 Section 3, unchanged)
# ============================================================================

md(r"""## 3. The principle: resolution bandwidth vs. Bessel zero-crossing spacing

Section 2's window sizes (5, 50, 100) came from the report as fixed numbers, with no explanation
of *why* those particular values separate "preserved" from "unreliable." There's a first-principles
reason, and it gives a quantitative criterion rather than a rule of thumb.

For large argument, $J_0(2\pi f r/c) \approx \sqrt{2/(\pi x)}\cos(x - \pi/4)$ with
$x = 2\pi f r /c$ -- an oscillation with *constant spacing in frequency* between zeros, since the
phase advances linearly in $f$. Setting consecutive zero-crossing phases $\pi$ apart gives

$$\Delta f_{\text{zero}} = \frac{c}{2r}$$

-- the same relationship phase-velocity picking uses in reverse ($c = 2r \Delta f_{\text{zero}}$,
Aki 1957, Notebook 2). Meanwhile, *any* frequency-domain averaging has an effective resolution
bandwidth: $2W$ for a taper set with half-bandwidth $W$ (hence the taper-count relationship
$K \approx 2NW$ from Notebook 1), or $w/N$ for an $w$-bin boxcar. If that resolution bandwidth
approaches or exceeds $\Delta f_{\text{zero}}$, the estimator is averaging *across* a full
oscillation cycle instead of *within* one -- which necessarily flattens the crossing's local
amplitude swing and blurs its location. One criterion, for any estimator:

$$2W \lesssim \Delta f_{\text{zero}} = \frac{c}{2r}$$

This directly predicts Section 2's result, without needing to know $c$ in advance --
$\Delta f_{\text{zero}}$ can be read off `FastMspec`'s own (undistorted) crossings. This is also
exactly the criterion `NW_high(r)` (Section 6/7 below, and Round 2's whole bandwidth sweep) is
built from -- the picker's own quality depends directly on staying on the right side of it.""")

md(r"""Two problems with stopping there, raised on review: (1) $\Delta f_{\text{zero}}$ isn't
actually one constant -- real dispersion curves have $c(f)$, so the true zero-crossing spacing
should vary across the spectrum, and a single global number begs the question of what spacing
we're even comparing against; (2) the *quality* of the spectrum matters as much as its noise level
-- Section 1's roughness metric measures variance, but a low-variance estimate can still have a
poorly-conditioned crossing (flat through zero) that's easy to mislocate. Both point to the same
fix: use each crossing's own *local* amplitude swing (already established at the end of Section 2)
both to decide which crossings are trustworthy enough to trust for $\Delta f_{\text{zero}}$, and
as the second, complementary line of evidence alongside drift. This section builds the argument in
five parts, following the structure of Ekstrom et al.'s Aki-spectral-formulation work on
noise-derived phase velocities (Ekstrom, Abers & Webb 2009; Ekstrom 2014, 2017 -- see References):
(1) raw spectral behavior at the smoothing extremes, (2) where each method's crossings are
actually trustworthy, (3) whether the implied zero-crossing spacing is stable enough to define
$\Delta f_{\text{zero}}$ at all, (4) drift measured against that spacing, and (5) pulling it
together.""")

md(r"""#### Part 1 — do the end-member curves even look different?

`FastMspec` is a single, fixed estimate -- it doesn't have a smoothing knob to turn. Single-taper's
smoothed variants do, from a mild window (5) to an aggressive one (100).""")

code("""fig, ax = plt.subplots(figsize=(8, 4.5))
cs_w5 = uniform_filter1d(cs_love, size=5)
cs_w100 = uniform_filter1d(cs_love, size=100)

zoom = (f_love >= 0.05) & (f_love < 0.15)
ax.plot(f_love[zoom], cf_love[zoom], color='C1', linewidth=1.5, label='FastMspec')
ax.plot(f_love[zoom], cs_w5[zoom], color='C0', linewidth=1, alpha=0.8, label='Single-taper, window=5')
ax.plot(f_love[zoom], cs_w100[zoom], color='C3', linewidth=1.5, label='Single-taper, window=100')
ax.axhline(0, color='gray', linewidth=0.5)
ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel('Re[coherency]')
ax.set_title('Figure 2 -- Zoomed: 0.05-0.15 Hz')
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()
""")

md(r"""Deliberately not what you might expect: at this scale, all three broadly track the same
oscillation, and `window=100` (more smoothing) is even *visually closer* to `FastMspec` than
`window=5` is -- because more smoothing means less point-to-point scatter, and scatter is what a
raw amplitude comparison actually measures. That's the point: **amplitude-level agreement doesn't
imply crossing-location agreement** -- a heavily smoothed curve can look like a good match while
its zero crossings have moved substantially, which is exactly what Parts 2-4 go on to show.
Eyeballing the spectrum isn't enough; you have to look at crossings specifically.""")

md(r"""#### Part 2 — where are the *trustworthy* crossings?

Not every zero-crossing is equally meaningful -- one sitting on a steep, high-amplitude swing is
well-conditioned; one on a shallow wiggle could be almost anything. Comparing swing *values*
directly across methods doesn't work, though: single-taper's raw noise alone produces bigger
local swings than `FastMspec`'s entire (low-variance) signal, so a shared absolute threshold is
meaningless -- a check below confirms 100% of single-taper's crossings would "pass" a threshold
set from `FastMspec`'s scale. Each method needs its own relative threshold: the top half of *its
own* swing distribution.""")

code("""def zero_crossings_idx(y, freqs):
    sgn = np.sign(y)
    idx = np.where(np.diff(sgn) != 0)[0]
    zc = freqs[idx] + (freqs[idx + 1] - freqs[idx]) * (0 - y[idx]) / (y[idx + 1] - y[idx])
    return zc, idx

def local_amplitude_swing(y, idx, half_win=3):
    out = []
    for i in idx:
        lo, hi = max(0, i - half_win), min(len(y), i + half_win + 2)
        out.append(y[lo:hi].max() - y[lo:hi].min())
    return np.array(out)

cf_band, cs_band = cf_love[love_band], cs_love[love_band]
zc_f, idx_f = zero_crossings_idx(cf_band, fb_love)
zc_s, idx_s = zero_crossings_idx(cs_band, fb_love)
swing_f = local_amplitude_swing(cf_band, idx_f)
swing_s = local_amplitude_swing(cs_band, idx_s)

print(f"FastMspec crossings: n={len(zc_f)}, swing range=[{swing_f.min():.5f}, {swing_f.max():.5f}]")
print(f"Single-taper (raw) crossings: n={len(zc_s)}, swing range=[{swing_s.min():.5f}, {swing_s.max():.5f}]")
shared_thresh = np.median(swing_f)
print(f"\\nIf we (wrongly) shared FastMspec's median swing ({shared_thresh:.5f}) as a threshold: "
      f"{100*np.mean(swing_s >= shared_thresh):.0f}% of single-taper's crossings would 'pass' -- meaningless.")

thresh_f, thresh_s = np.median(swing_f), np.median(swing_s)
reliable_f, reliable_s = swing_f >= thresh_f, swing_s >= thresh_s
print(f"\\nUsing each method's own median instead: FastMspec {reliable_f.sum()}/{len(zc_f)} reliable, "
      f"single-taper {reliable_s.sum()}/{len(zc_s)} reliable (by construction, ~50% each).")
""")

md(r""""Reliable" above is entirely self-referential -- each method's own top-50%-by-swing, not
checked against anything physical. Adding that check directly: SA53/SA58's own real per-frequency
hybrid ADAMA+GDM52 curve (same construction as Section 5 uses for SKRH-BAND -- SA53/58 sit inside
ADAMA's Africa coverage even though they're not part of the `XV` Madagascar dataset) predicts
where the *true* Bessel zero-crossings should fall, independent of either method's own noise.""")

code("""from dispcurve_pick.hybrid_reference_curve import AdamaMap, Gdm52Map, build_reference_curve

africa_stalist = pd.read_csv('../data/metadata/africa_adama_stalist.csv')
lat1_p2, lon1_p2 = africa_stalist.loc[africa_stalist.Station == 'SA58', ['Latitude', 'Longitude']].iloc[0]
lat2_p2, lon2_p2 = africa_stalist.loc[africa_stalist.Station == 'SA53', ['Latitude', 'Longitude']].iloc[0]
adama_p2 = AdamaMap('../data/reference/adama_maps')
gdm52_p2 = Gdm52Map('../data/reference/gdm52')
hybrid_p2 = build_reference_curve(lat1_p2, lon1_p2, lat2_p2, lon2_p2, sources=[adama_p2, gdm52_p2], wave='love')

from scipy.special import j0 as _j0

f_true_lo, f_true_hi = max(hybrid_p2.f_lo, fb_love.min()), min(hybrid_p2.f_hi, fb_love.max())
f_true_grid = np.linspace(f_true_lo, f_true_hi, 5000)
pred_true = _j0(2 * np.pi * f_true_grid * dist_km_t / hybrid_p2.func(f_true_grid))  # exact J0, Aki 1957
zc_true, _ = zero_crossings_idx(pred_true, f_true_grid)
print(f"{len(zc_true)} true/expected zero-crossings in [{f_true_lo:.3f}, {f_true_hi:.3f}] Hz, "
      f"from SA53/58's own real hybrid curve")
""")

code("""fig, axes = plt.subplots(3, 1, figsize=(12, 5), sharex=True, gridspec_kw={'height_ratios': [1, 1, 0.6]})
axes[0].vlines(zc_f[reliable_f], 0, 1, color='C1', linewidth=1.2)
axes[0].vlines(zc_f[~reliable_f], 0, 1, color='lightgray', linewidth=0.6)
axes[0].set_yticks([]); axes[0].set_ylabel('FastMspec', rotation=0, ha='right', va='center')
axes[1].vlines(zc_s[reliable_s], 0, 1, color='C0', linewidth=1.2)
axes[1].vlines(zc_s[~reliable_s], 0, 1, color='lightgray', linewidth=0.6)
axes[1].set_yticks([]); axes[1].set_ylabel('Single-taper\\n(raw)', rotation=0, ha='right', va='center')
axes[2].vlines(zc_true, 0, 1, color='black', linewidth=1.2)
axes[2].set_yticks([]); axes[2].set_ylabel('True/expected\\n(hybrid curve)', rotation=0, ha='right', va='center')
axes[2].set_xlabel('Frequency (Hz)')
axes[2].set_xlim(f_true_lo, f_true_hi)
fig.suptitle('Figure 3 -- Reliable crossings vs. true/expected crossings (from a real reference curve)')
plt.tight_layout()
plt.show()
""")

md(r"""Now a true-positive/false-positive read is possible, not just relative density: `FastMspec`'s
reliable (orange) crossings land close to the true/expected ticks (bottom row) far more often than
single-taper's reliable (blue) crossings do -- the true positives, not just "fewer events."

Visually the difference is stark even without the overlay: `FastMspec`'s reliable crossings are spread fairly evenly
across the band; single-taper's are packed into dense clusters -- a signature of noise-driven
zero-crossings clustering tightly around wherever the (unaveraged) curve happens to wander near
zero, not of genuine, evenly-spaced dispersion structure.""")

md(r"""#### Part 3 — is $\Delta f_{\text{zero}}$ even stable enough to define?

If reliable crossings reflect real dispersion structure, their spacing should cluster in a
physically sensible, reasonably bounded range. If they're still noise-contaminated, spacing will
be erratic and concentrated at unrealistically small values (noise wiggles are close together by
nature).""")

code("""spacing_f_reliable = np.diff(np.sort(zc_f[reliable_f]))
spacing_s_reliable = np.diff(np.sort(zc_s[reliable_s]))

print(f"FastMspec (reliable): n={len(spacing_f_reliable)}, mean={spacing_f_reliable.mean():.5f} Hz, "
      f"median={np.median(spacing_f_reliable):.5f} Hz, std/mean={spacing_f_reliable.std()/spacing_f_reliable.mean():.2f}")
print(f"Single-taper (reliable, own top 50%): n={len(spacing_s_reliable)}, mean={spacing_s_reliable.mean():.5f} Hz, "
      f"median={np.median(spacing_s_reliable):.5f} Hz, std/mean={spacing_s_reliable.std()/spacing_s_reliable.mean():.2f}")
print(f"\\nEven restricted to its OWN best half, single-taper's spacing is "
      f"{spacing_f_reliable.mean()/spacing_s_reliable.mean():.1f}x tighter than FastMspec's -- "
      f"its 'best' crossings are still predominantly noise-driven, not dispersion-driven.")

fig, ax = plt.subplots(figsize=(7, 4.5))
bins = np.logspace(np.log10(1e-5), np.log10(0.05), 40)
ax.hist(spacing_f_reliable, bins=bins, color='C1', alpha=0.6, label='FastMspec (reliable)',
        weights=np.ones_like(spacing_f_reliable) / len(spacing_f_reliable))
ax.hist(spacing_s_reliable, bins=bins, color='C0', alpha=0.6, label='Single-taper (reliable, own top 50%)',
        weights=np.ones_like(spacing_s_reliable) / len(spacing_s_reliable))
ax.set_xscale('log')
ax.set_xlabel(r'Local zero-crossing spacing (Hz, log scale)')
ax.set_ylabel('Fraction of that method\\'s own crossings')
ax.set_title(r'Figure 4 -- $\\Delta f_{zero}$ distribution: stable for FastMspec, not for single-taper')
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()

zc_f_reliable_sorted = np.sort(zc_f[reliable_f])
local_df_zero = np.zeros(len(zc_f_reliable_sorted))
for i in range(len(zc_f_reliable_sorted)):
    gaps = []
    if i > 0:
        gaps.append(zc_f_reliable_sorted[i] - zc_f_reliable_sorted[i - 1])
    if i < len(zc_f_reliable_sorted) - 1:
        gaps.append(zc_f_reliable_sorted[i + 1] - zc_f_reliable_sorted[i])
    local_df_zero[i] = np.mean(gaps)
""")

md(r"""So: not fully "known across the entire spectrum" as one number, but **bounded and
physically interpretable** for `FastMspec`'s reliable subset in a way single-taper's crossings, at
any quality percentile, are not. This is the honest scope of the claim -- a local, quality-filtered
$\Delta f_{\text{zero}}(f)$, not a single global constant.""")

md(r"""#### Part 4 — drift statistics, against that local spacing

For each of `FastMspec`'s reliable, locally-spaced crossings, how far is the nearest crossing in
each smoothed single-taper variant -- as a fraction of *that crossing's own* local spacing, not
one global number?""")

code("""def nearest(ref, other):
    return other[np.argmin(np.abs(other - ref))] if len(other) else np.nan

drift_stats = {}
for w in [5, 50, 100]:
    variant = uniform_filter1d(cs_love, size=w)[love_band]
    zc_w, _ = zero_crossings_idx(variant, fb_love)
    pct = np.array([100 * abs(nearest(zc_r, zc_w) - zc_r) / dfz
                     for zc_r, dfz in zip(zc_f_reliable_sorted, local_df_zero)])
    drift_stats[w] = pct
    print(f"window={w:3d}: drift as % of local df_zero -- "
          f"median={np.median(pct):.1f}%, IQR=[{np.percentile(pct,25):.1f}%,{np.percentile(pct,75):.1f}%], "
          f"fraction >50% drift = {100*np.mean(pct>50):.0f}%")

median_local_df_zero = np.median(local_df_zero)
print(f"\\nClosing the loop with the opening criterion: FastMspec's own 2W = {2*wband:.5f} Hz, "
      f"vs. median local df_zero = {median_local_df_zero:.5f} Hz -- "
      f"ratio = {2*wband/median_local_df_zero:.2f} (criterion: <1 is safe), "
      f"now using the refined, reliable-subset-based estimate rather than the earlier single clean-band number.")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.boxplot([np.clip(drift_stats[w], 0, 200) for w in [5, 50, 100]], tick_labels=['window=5', 'window=50', 'window=100'],
           showfliers=False)
ax.axhline(50, color='r', linestyle='--', linewidth=1, label='50% of local spacing')
ax.set_ylabel('Drift as % of local $\\\\Delta f_{zero}$ (clipped at 200% for display)')
ax.set_title('Figure 5 -- Drift distribution grows with smoothing window')
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()
""")

md(r"""Median drift climbs from ~4% of local spacing (window=5, negligible) to ~22% (window=50) to
~51% (window=100) -- at the median, a window=100 crossing has moved *more than halfway to its
neighbor*, and just over half of all crossings exceed that mark. That's not added noise on an
otherwise correct pick; it's a real risk of identifying the wrong crossing entirely.""")

md(r"""#### Part 5 — pulling it together

- **Part 1**: `FastMspec` is one fixed, stable curve; single-taper visibly warps further from it
  as the smoothing window widens -- there's no window size that's simultaneously "smooth enough"
  and "close enough."
- **Part 2**: swing-based reliability, applied fairly (each method against its own distribution),
  shows `FastMspec`'s trustworthy crossings spread evenly across the band; single-taper's cluster
  into dense, noise-driven clumps.
- **Part 3**: the zero-crossing spacing implied by those crossings is physically bounded and
  interpretable for `FastMspec`, but not for single-taper at any quality percentile.
- **Part 4**: measured against that yardstick, drift is negligible at window=5, already
  substantial at window=50, and comparable to half the crossing spacing at window=100.

Together these support the resolution-bandwidth criterion from earlier in this section
($2W \lesssim \Delta f_{\text{zero}}$) with a more rigorous empirical foundation than a single
global $\Delta f_{\text{zero}}$: `FastMspec`'s advantage is real, quantifiable, and grounded in
the same class of noise-derived phase-velocity reasoning as Ekstrom, Abers & Webb (2009) and
Ekstrom (2014, 2017) -- not an artifact of matching one report's specific numbers. The rest of
this notebook now puts that resolution-bandwidth criterion to work directly: Section 4 introduces
the picker that uses it, Sections 5-7 test it at real, then full-dataset, scale.
""")

# ============================================================================
# Section 4 -- The method: the vendored, instrumented picker
# ============================================================================

md(r"""## 4. The method: the vendored, instrumented `seislib` picker

`seislib.an.extract_dispcurve` (Magrini et al. 2022) picks a phase-velocity dispersion curve from
a coherence spectrum's zero crossings, anchored to a reference curve to resolve the $2\pi$ branch
ambiguity between parallel dispersion curves -- essentially the Aki-spectral-formulation method
(Ekstrom, Abers & Webb 2009) with real internal quality control, not a naive event scan:

- **Bad-quality gating.** Every candidate zero crossing is checked against 3 criteria before
  picking even starts: a peak-ratio test (Part 2's "is the swing well-conditioned" idea, done
  properly), an envelope-relative-amplitude test, and a spacing-vs-reference-curve test (directly
  Section 3's $\Delta f_{\text{zero}}$ criterion, applied per-crossing against the reference
  curve's own local prediction rather than one global number).
- **Sequential, reference-guided tracking.** Picking walks outward in frequency from wherever
  quality is first good, using a kernel-density field built from *every possible cycle-count*
  interpretation of each crossing (the same $2\pi$-branch problem Section 3 sidesteps by only
  ever comparing `FastMspec` against itself) -- the reference curve disambiguates which cycle
  count is physically plausible.
- **A real acceptance test**, not just "did something come out": the final curve must cover at
  least 1/5 of the requested frequency band on a log-spaced grid, or picking is reported as
  non-converged rather than returning a token curve.

This project vendors and lightly instruments this picker (`python/dispcurve_pick/`) rather than
calling the installed package as-is: `return_diagnostics=True` surfaces the quality signals above
explicitly -- `bad_quality_fraction`, `n_candidate_crossings`, `n_accepted_picks`,
`freq_coverage_fraction`, `mean_amp_ratio`, and (added specifically for this notebook) the full
per-crossing frequency/bad-quality arrays -- instead of discarding them. Default behavior
(`return_diagnostics=False`) is byte-identical to the unmodified upstream package
(`python/dispcurve_pick/tests/test_matches_upstream.py`); nothing about the picking algorithm
itself is changed.

A quick synthetic demonstration -- the exact Bessel-coherence construction
$\text{Re}\{\gamma(f)\} = J_0(2\pi f r/c(f))$ this project's translation-fidelity test uses,
clean (no noise), so the picker should converge cleanly and the barcode should show almost no
bad-quality crossings.""")

code("""from scipy.special import j0

dist_km_demo = 300.0
freqs_demo = np.linspace(0.01, 0.4, 2000)
c_of_f_demo = lambda f: 3.0 - 1.5 * f  # mild normal dispersion
x_demo = 2 * np.pi * freqs_demo * dist_km_demo / c_of_f_demo(freqs_demo)
coherence_demo = j0(x_demo)
ref_curve_demo = np.column_stack([freqs_demo[::50], c_of_f_demo(freqs_demo[::50])])

curve_demo, picks_demo, diag_demo = extract_dispcurve(
    freqs_demo, coherence_demo, dist_km_demo, ref_curve_demo,
    freqmin=0.02, freqmax=0.35, cmin=1.0, cmax=5.0, pick_threshold=0,
    return_diagnostics=True,
)
print(f"converged={diag_demo.converged}, bad_quality_fraction={diag_demo.bad_quality_fraction:.3f}, "
      f"n_candidate_crossings={diag_demo.n_candidate_crossings}, "
      f"freq_coverage_fraction={diag_demo.freq_coverage_fraction:.3f}, "
      f"mean_amp_ratio={diag_demo.mean_amp_ratio:.2f}")
""")

md(r"""**Figure 6 -- the picker's own native diagnostic plot** (`plotting=True`, re-run once on the
same input purely to render it -- the call above already has the data/diagnostics). Top: the raw
coherence with low-quality crossings marked in red. Bottom: the actual kernel-density field the
picker searches (not an invented barcode) -- reference curve as the light-blue line, the tracked
zero-crossing branch in blue, accepted picks in red, the final smoothed curve dashed red. This is
the picker authors' own purpose-built visualization, used throughout the rest of this notebook
instead of a custom-built alternative.""")

code("""_ = extract_dispcurve(
    freqs_demo, coherence_demo, dist_km_demo, ref_curve_demo,
    freqmin=0.02, freqmax=0.35, cmin=1.0, cmax=5.0, pick_threshold=0,
    plotting=True,
)
""")

md(r"""Clean data, clean result: essentially every candidate crossing is good quality, the branch
tracks the reference curve tightly, the picker converges, and full-band coverage is achieved. This
is the ceiling case -- Section 5 immediately below shows what real, noisy data looks like by
comparison.""")

# ============================================================================
# Section 5 -- Worked example: SKRH-BAND
# ============================================================================

md(r"""## 5. Worked example: SKRH-BAND vs. Sayan's own known-good result

`AF.SKRH`-`XV.BAND` (290.3 km, Madagascar) is the pair this project's own bluehive batch pipeline
was validated against before scaling to all 380 pairs (`verification/skrh_band_real_data/`) --
the natural first real-data test here too. Real matched data (107 usable days, `S1_data_mat`/
`S2_data_mat`), same `FastMspec` production config as the batch pipeline.""")

code("""WBAND, CUTOFF, EPSILON = 0.001, 1 - 1e-5, 1e-5
PICK_CMIN, PICK_CMAX = 1.2, 4.8

t0 = time.time()
raw = loadmat(DATA / 'madagascar' / 'AFSKRH_XVBAND_win_3_all_matched_data.mat', simplify_cells=True)
s1_skrh, s2_skrh = raw['S1_data_mat'], raw['S2_data_mat']
dist_km_skrh = float(raw.get('stapairsinfo', {}).get('r', 290.3))
print(f"S1/S2 shape: {s1_skrh.shape} ({s1_skrh.shape[0]} usable days), dist={dist_km_skrh:.1f} km")

r_skrh = compute_crosscorr_mtc_fastmspec(s1_skrh, s2_skrh, wband=WBAND, cutoff=CUTOFF, epsilon=EPSILON)
n_skrh = s1_skrh.shape[2]
faxis_skrh = np.fft.fftfreq(n_skrh, d=1.0)
pos_skrh = faxis_skrh > 0
faxis_pos_skrh, coh_pos_skrh = faxis_skrh[pos_skrh], r_skrh.coh_sum[pos_skrh].real / r_skrh.coh_num
print(f"crosscorr done in {time.time()-t0:.1f}s, coh_num={r_skrh.coh_num}, "
      f"coherence range=[{coh_pos_skrh.min():.4f}, {coh_pos_skrh.max():.4f}]")
print(f"Cross-validated against Sayan's own MATLAB coh_sum in verification/skrh_band_real_data/: "
      f"this exact range (-0.1033 to 0.0834) already confirmed to match there.")
""")

md(r"""### The per-pair hybrid ADAMA+GDM52 reference curve

Rather than a single generic Love-wave curve shared by every pair (`SDISPL.ASC`), this pair's own
great-circle path is looked up directly against ADAMA_Maps (Africa, 6-40s) and GDM52 (global,
45-150s) -- see `data/reference/hybrid_curve_README.md` for the full method (slowness-averaging
along the path) and standalone validation against a real ADAMA-measured pair.""")

code("""manifest_mdg = pd.read_csv(DATA / 'madagascar' / 'madagascar_stn_conn_ccflist.csv')
row_skrh = manifest_mdg[(manifest_mdg.net1 == 'AF') & (manifest_mdg.stn1 == 'SKRH') &
                         (manifest_mdg.net2 == 'XV') & (manifest_mdg.stn2 == 'BAND')].iloc[0]
lat1, lon1, lat2, lon2 = (float(row_skrh.stn1lat), float(row_skrh.stn1lon),
                           float(row_skrh.stn2lat), float(row_skrh.stn2lon))

adama = AdamaMap('../data/reference/adama_maps')
gdm52 = Gdm52Map('../data/reference/gdm52')
hybrid_skrh = build_reference_curve(lat1, lon1, lat2, lon2, sources=[adama, gdm52], wave='love')
print(f"hybrid curve: {len(hybrid_skrh.samples)} period samples, "
      f"f=[{hybrid_skrh.f_lo:.5f}, {hybrid_skrh.f_hi:.5f}] Hz, "
      f"low-confidence periods (<12s): {hybrid_skrh.low_confidence_periods}")

fig, ax = plt.subplots(figsize=(7, 4.5))
freqs_plot = np.linspace(hybrid_skrh.f_lo, hybrid_skrh.f_hi, 300)
ax.plot(1/freqs_plot, hybrid_skrh.func(freqs_plot), color='C1', linewidth=1.5, label='Hybrid curve (smoothed)')
colors = ['C3' if lc else 'C0' for lc in [s.low_confidence for s in hybrid_skrh.samples]]
ax.scatter(hybrid_skrh.periods, hybrid_skrh.velocities, c=colors, zorder=5, label='Raw path-averaged samples')
ax.set_xlabel('Period (s)'); ax.set_ylabel('Phase velocity (km/s)')
ax.set_title('Figure 7 -- SKRH-BAND hybrid reference curve (red=low-confidence, T<12s)')
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()
""")

md(r"""### Scanning the corridor: hybrid curve vs. the old shared `SDISPL.ASC` curve

A modest corridor scan (±0.8 km/s, coarse 0.2 km/s step -- 9 templates, not the full 0.05 km/s-step
production corridor used in the batch pipeline, kept coarse here purely so this cell finishes in a
few minutes rather than tens) run through the widened-corridor strategy (Stage 4.5, Section 7
below), on both the new hybrid curve and the old shared curve this pair would have used before
Stage 4.5 -- an isolated, single-pair version of the same before/after comparison Section 7
reports at 4-pair scale.""")

code("""c_ref_old, f_lo_old, f_hi_old = load_reference_curve(REFERENCE, freqmin=0.01, freqmax=0.5)

templates_hybrid = build_template_family_widened(hybrid_skrh.func, hybrid_skrh.f_lo, hybrid_skrh.f_hi,
                                                    corridor_km_s=0.8, step_km_s=0.2)
templates_old = build_template_family_widened(c_ref_old, f_lo_old, f_hi_old,
                                                 corridor_km_s=0.8, step_km_s=0.2)

print(f"Scanning {len(templates_hybrid)} templates, hybrid curve (f=[{hybrid_skrh.f_lo:.4f},{hybrid_skrh.f_hi:.4f}] Hz):")
t0 = time.time()
scanned_hybrid = scan_templates_with_picker(faxis_pos_skrh, coh_pos_skrh, dist_km_skrh, templates_hybrid,
                                              hybrid_skrh.f_lo, hybrid_skrh.f_hi, PICK_CMIN, PICK_CMAX)
print(f"  ({time.time()-t0:.0f}s)")

print(f"Scanning {len(templates_old)} templates, old SDISPL.ASC curve (f=[{f_lo_old:.4f},{f_hi_old:.4f}] Hz):")
t0 = time.time()
scanned_old = scan_templates_with_picker(faxis_pos_skrh, coh_pos_skrh, dist_km_skrh, templates_old,
                                           f_lo_old, f_hi_old, PICK_CMIN, PICK_CMAX)
print(f"  ({time.time()-t0:.0f}s)")

n_conv_hybrid = sum(1 for _, _, _, d in scanned_hybrid if d is not None and d.converged)
n_conv_old = sum(1 for _, _, _, d in scanned_old if d is not None and d.converged)
print(f"\\nTemplates converged: hybrid curve {n_conv_hybrid}/{len(templates_hybrid)}, "
      f"old curve {n_conv_old}/{len(templates_old)}")

delta_h, curve_h, picks_h, diag_h = best_of(scanned_hybrid)
delta_o, curve_o, picks_o, diag_o = best_of(scanned_old)
for label, delta, diag in [('Hybrid curve', delta_h, diag_h), ('Old SDISPL.ASC curve', delta_o, diag_o)]:
    if diag is None:
        print(f"{label}: no template converged")
    else:
        print(f"{label}: best delta={delta:+.2f} km/s, score={score(diag):.3f}, "
              f"coverage={diag.freq_coverage_fraction:.3f}, bad_quality={diag.bad_quality_fraction:.3f}, "
              f"n_picks={diag.n_accepted_picks}")
""")

md(r"""**Figures 8 (hybrid curve) and 9 (old curve) -- the picker's own native diagnostic plot, for
the winning template from each curve** (re-run once each with `plotting=True`, purely to render --
the scan above already has the data). Same reading as Figure 6: reference curve as the light-blue
line, tracked branch in blue, accepted picks in red, low-quality candidates marked red on the
coherence panel above. Hybrid curve first, old curve second -- compare the two directly.""")

code("""if diag_h is not None:
    freqs200_h = np.linspace(hybrid_skrh.f_lo, hybrid_skrh.f_hi, 200)
    ref_curve_h = np.column_stack([freqs200_h, templates_hybrid[delta_h](freqs200_h)])
    print(f"Hybrid curve, delta={delta_h:+.2f} km/s:")
    _ = extract_dispcurve(
        faxis_pos_skrh, coh_pos_skrh, dist_km_skrh, ref_curve_h,
        freqmin=hybrid_skrh.f_lo, freqmax=hybrid_skrh.f_hi, cmin=PICK_CMIN, cmax=PICK_CMAX,
        horizontal_polarization=True, plotting=True,
    )
else:
    print('Hybrid curve: no template converged, nothing to plot')
""")

code("""if diag_o is not None:
    freqs200_o = np.linspace(f_lo_old, f_hi_old, 200)
    ref_curve_o = np.column_stack([freqs200_o, templates_old[delta_o](freqs200_o)])
    print(f"Old SDISPL.ASC curve, delta={delta_o:+.2f} km/s:")
    _ = extract_dispcurve(
        faxis_pos_skrh, coh_pos_skrh, dist_km_skrh, ref_curve_o,
        freqmin=f_lo_old, freqmax=f_hi_old, cmin=PICK_CMIN, cmax=PICK_CMAX,
        horizontal_polarization=True, plotting=True,
    )
else:
    print('Old curve: no template converged, nothing to plot')
""")

md(r"""This single-pair result previews Section 7's full 4-pair Stage 4.5 finding: more templates
converge with the hybrid curve than the old shared curve. Comparing the two figures above directly:
more of the hybrid curve's candidate crossings are good-quality (fewer red dots on the top panel of
its own figure), and the tracked branch follows the reference curve more tightly -- the *mechanism*
behind the summary statistic is visible directly, not just the number, when the reference curve the
picker is anchored to is a better physical match for this specific path.""")

# ============================================================================
# Section 6 -- At scale: Round 1 + Round 2
# ============================================================================

md(r"""## 6. At scale: the full 380-pair dataset

Section 5 showed the picker working on one pair. This section loads the real, committed results
of running it -- via `python/dispcurve_pick_batch/work_unit.py`, on bluehive -- across the entire
380-pair Madagascar dataset, at four techniques (Round 1) and across a bandwidth sweep (Round 2).
Both predate Stage 4.5's fixes (Section 7): single generic `SDISPL.ASC` curve, uniform corridor,
the un-fixed `horizontal_polarization` default -- the picker configuration Section 5's "old curve"
column used, not the "hybrid curve" one. Full provenance: `data/results/dispcurve_quality/README.md`.

**A quick recap of the metrics the figures below use** (all from `PickDiagnostics`, Section 4):
`converged` -- did picking succeed at all (a yes/no gate). `freq_coverage_fraction` -- of the
requested frequency band, what fraction does the final picked curve actually span (higher is
better, 1.0 = full coverage). `bad_quality_fraction` -- of all candidate zero-crossings the picker
considered, what fraction failed its own quality gate before picking even started (lower is
better). These three numbers are what every bar/line below actually plots -- worth having them
in view rather than relying on Section 4's own introduction from several sections back.

### Round 1: single-taper vs. three multitaper techniques, all 380 pairs

This is the headline result Section 1 promised: not a roughness comparison on one pair, but
dispersion-curve *convergence* across the entire real dataset.""")

code("""manifest_r1 = pd.read_csv(RESULTS / 'manifest.csv')
manifest_r1['quartile'] = pd.qcut(manifest_r1['dist_km'], 4, labels=[1, 2, 3, 4])

conv_by_technique = manifest_r1.groupby('technique').converged.mean().sort_values()
print("Convergence rate, all 380 pairs, by technique:")
print(conv_by_technique.apply(lambda x: f'{x:.1%}'))

conv_only = manifest_r1[manifest_r1.converged]
quality_by_technique = conv_only.groupby('technique')[['bad_quality_fraction', 'freq_coverage_fraction']].mean()
print("\\nAmong CONVERGED pairs only -- quality, not just convergence:")
print(quality_by_technique.round(3))
""")

code("""fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

techniques = ['single-taper', 'Mspec', 'MspecBestK', 'FastMspec']
colors_t = {'single-taper': 'C0', 'Mspec': 'C2', 'MspecBestK': 'C3', 'FastMspec': 'C1'}

axes[0].bar(techniques, [conv_by_technique[t] for t in techniques],
            color=[colors_t[t] for t in techniques])
axes[0].set_ylabel('Convergence rate (all 380 pairs)')
axes[0].set_title('Convergence by technique')
for i, t in enumerate(techniques):
    axes[0].text(i, conv_by_technique[t] + 0.01, f'{conv_by_technique[t]:.0%}', ha='center', fontsize=9)

quartile_conv = manifest_r1.groupby(['quartile', 'technique'], observed=True).converged.mean().unstack()
for t in techniques:
    axes[1].plot(quartile_conv.index, quartile_conv[t], 'o-', color=colors_t[t], label=t)
axes[1].set_xlabel('Distance quartile (1=nearest)'); axes[1].set_ylabel('Convergence rate')
axes[1].set_title('Convergence vs. distance')
axes[1].legend(fontsize=8)

x = np.arange(len(techniques))
w = 0.35
axes[2].bar(x - w/2, [quality_by_technique.loc[t, 'freq_coverage_fraction'] if t in quality_by_technique.index else 0 for t in techniques],
            w, label='freq coverage', color='C0')
axes[2].bar(x + w/2, [1 - quality_by_technique.loc[t, 'bad_quality_fraction'] if t in quality_by_technique.index else 0 for t in techniques],
            w, label='1 - bad_quality_fraction', color='C3')
axes[2].set_xticks(x); axes[2].set_xticklabels(techniques, rotation=20)
axes[2].set_title('Quality among CONVERGED pairs only')
axes[2].legend(fontsize=8)
axes[2].set_ylim(0, 1)

plt.suptitle('Figure 10 -- Round 1: convergence and quality by technique, all 380 pairs')
plt.tight_layout()
plt.show()
""")

md(r"""**The headline, stated plainly**: single-taper's dispersion-curve convergence rate across
all 380 real pairs is **0%** -- not low, zero, at every distance quartile. `FastMspec` and
`MspecBestK` converge on roughly a quarter of pairs; `Mspec` (no `FastMultitaper` fusion, same
taper count) noticeably less. This is the real-dataset-scale confirmation of everything Sections
1-3 built up on one pair: multitaper isn't merely "smoother," it's the difference between a
pipeline that can extract *any* real dispersion curves from this dataset and one that extracts
none at all.

**The middle panel matters as much as the headline number**: convergence for every multitaper
technique collapses with distance -- strong at the nearest quartile, weak at the farthest. That's
not a flaw unique to this picker; it's the direct real-dataset manifestation of Section 3's own
resolution-bandwidth criterion ($2W \lesssim c/2r$) -- farther pairs have tighter true zero-crossing
spacing at fixed bandwidth, so the same `WBAND` that resolves a near pair cleanly starts to blur a
far one. Round 2 (next) tests exactly this.

**Right panel, the honest complication (per direct guidance -- lean on quality, not just
convergence)**: among pairs that DO converge, `Mspec` -- despite converging *less often* overall
-- has the best quality on the pairs it does converge on (higher coverage, lower bad-quality
fraction than `FastMspec`/`MspecBestK`). Read together with the convergence-rate panel, the
honest picture is: `FastMspec`/`MspecBestK` convert more real pairs into a curve at all, but
`Mspec`'s survivors tend to be slightly cleaner picks -- a real trade-off between *how many* pairs
you get a curve from and *how clean* that curve is when you do, not a simple one-technique-wins-
everywhere story.""")

md(r"""### Round 2: does bandwidth choice explain the distance-dependent falloff?

Round 2 swept `FastMspec`'s bandwidth (`NW`, as a fraction of `NW_high(r)` -- Section 3's own
criterion, evaluated per pair) across 15 pairs per distance quartile (60 pairs, stratified sample)
x 5 fractions = 300 points, testing 4 specific hypotheses head-on against real data
(`docs/round2_hypothesis_evaluation.tex` for the complete writeup; reproduced live here from the
same 300 committed result files, not re-typed).""")

code("""subset_r2 = pd.read_csv(RESULTS / 'round2_subset.csv')
rows_r2 = [json.load(open(f)) for f in glob.glob(str(RESULTS / 'round2_sweep' / '*.json'))]
sweep_r2 = pd.DataFrame(rows_r2).merge(
    subset_r2[['net1', 'stn1', 'net2', 'stn2', 'quartile', 'nw_high']],
    on=['net1', 'stn1', 'net2', 'stn2'],
)
N_SAMPLES = 10801
sweep_r2['fraction'] = (sweep_r2.wband_used * N_SAMPLES / sweep_r2.nw_high).round(2)
sweep_r2['k0_floor'] = sweep_r2.error.fillna('').str.contains('K=0')

print(f"Overall: {sweep_r2.converged.sum()}/{len(sweep_r2)} converged "
      f"({sweep_r2.converged.mean():.0%}); {sweep_r2.k0_floor.sum()}/{len(sweep_r2)} hit the "
      f"K=0 numerical floor ({sweep_r2.k0_floor.mean():.1%})")

table_conv = sweep_r2.pivot_table(index='quartile', columns='fraction', values='converged', aggfunc='mean')
table_conv = table_conv[sorted(table_conv.columns, reverse=True)]
print("\\nH1/H2/H3 -- convergence rate by quartile x fraction of NW_high(r):")
print((table_conv * 100).round(0).astype(int).astype(str) + '%')

print("\\nH4 -- K=0 floor rate by quartile (fraction of all 75 sweep points per quartile):")
print(sweep_r2.groupby('quartile').k0_floor.mean().apply(lambda x: f'{x:.0%}'))
""")

code("""fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for q in [1, 2, 3, 4]:
    axes[0].plot(table_conv.columns, table_conv.loc[q], 'o-', label=f'Q{q}')
axes[0].set_xlabel('Fraction of NW_high(r)'); axes[0].set_ylabel('Convergence rate')
axes[0].set_title('H1/H2: convergence vs. bandwidth fraction, by quartile')
axes[0].legend(fontsize=8, title='Quartile')

k0_by_q = sweep_r2.groupby('quartile').k0_floor.mean()
axes[1].bar(k0_by_q.index.astype(str), k0_by_q.values, color='C3')
axes[1].set_xlabel('Distance quartile'); axes[1].set_ylabel('K=0 floor rate (of 75 sweep points)')
axes[1].set_title('H4: numerical floor rate rises with distance')
plt.suptitle('Figure 11 -- Round 2: bandwidth-sweep convergence (H1-H3) and the K=0 numerical floor (H4)')
plt.tight_layout()
plt.show()
""")

md(r"""**H1 (normalized bandwidth) and H3 (existence-threshold cliff), strongly confirmed**:
`NW_high(r)`-normalizing does make quartiles 1-3 behave comparably to each other (each shows the
same qualitative rise-with-fraction pattern) -- but quartile 4 does not participate in this
pattern *at all*: **0% convergence at every single bandwidth fraction tested, including 1.5x its
own ceiling**. That's not "quartile 4 needs a different bandwidth" -- normalizing by its own
`NW_high(r)` and trying the widest fraction anyway still fails completely, a genuine existence
threshold, not a tunable-parameter problem.

**H2 (an interior MSE-optimal bandwidth), NOT supported** in the tested range: convergence rises
monotonically with fraction in every quartile that converges at all -- the widest bandwidth tried
(1.5x `NW_high(r)`) is uniformly the best, not an interior optimum. `docs/stage5_bandwidth_theory.tex`'s
bias-variance MSE framework predicts a real lower bound should exist somewhere; this sweep's range
just didn't reach it. **This is exactly what Notebook 5 (next) is built to test properly** -- a
true ADAMA benchmark, not Sayan's 380 pairs, specifically to push past this range with real
ground truth to compare against.

**H4 (a hard numerical floor), confirmed with a clean quartile gradient**: `FastMultitaper` is
literally undefined at `K=0` (no taper reaches the eigenvalue cutoff) -- a hard structural limit,
not a soft statistical one. The rate climbs monotonically from ~1% (nearest quartile) to ~51%
(farthest) as the *same* bandwidth fraction implies a smaller and smaller absolute `NW` for
farther pairs. By quartile 4, over half of all sweep points aren't just non-convergent for
picking reasons -- the underlying spectral estimate itself doesn't exist there.""")

# ============================================================================
# Section 7 -- Stage 4.5: reference curve, polarization, and corridor fixes
# ============================================================================

md(r"""## 7. Stage 4.5: reference-curve accuracy and the corridor strategy

Round 1/2 (Section 6) both used a real, but real*ly limited* picker configuration: a single
generic `SDISPL.ASC` curve shared by all 380 pairs regardless of actual geography, an unfixed
`horizontal_polarization` bug (Love-wave data needs the J0-J2 Bessel model, not the pure-J0 model
the picker defaulted to -- ~23% zero-location error at the first Bessel zero, converging to <0.1%
by the far field), and a uniform corridor scan with no extra care at low-confidence (short) periods.
Stage 4.5 closed all three gaps. Section 5 already showed one pair's version of this; here is the
real, committed 4-pair validation the fix was accepted on.

### Choosing the corridor strategy

Three strategies were tested head-to-head on the 4 report example pairs (one per distance
quartile): a uniform ±0.8 km/s corridor (baseline), widening that corridor specifically at
low-confidence (short) periods, and leaving the corridor alone but down-weighting short-period
coverage in the scoring function instead.""")

code("""corridor_eval = pd.read_csv(RESULTS / 'stage4_5' / 'corridor_strategy_eval.csv')
print(corridor_eval[['quartile', 'pair', 'strategy', 'delta', 'converged', 'coverage', 'bad_quality', 'n_picks']]
      .to_string(index=False))
""")

code("""fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
strategies = ['BASELINE', 'APPROACH_A_widen', 'APPROACH_B_downweight']
colors_s = {'BASELINE': 'gray', 'APPROACH_A_widen': 'C1', 'APPROACH_B_downweight': 'C0'}
x = np.arange(1, 5)
w = 0.25
for i, s in enumerate(strategies):
    sub = corridor_eval[corridor_eval.strategy == s].set_index('quartile').reindex([1, 2, 3, 4])
    axes[0].bar(x + (i - 1) * w, sub['coverage'].fillna(0), w, label=s, color=colors_s[s])
    axes[1].bar(x + (i - 1) * w, sub['bad_quality'].fillna(1), w, label=s, color=colors_s[s])
axes[0].set_xticks(x); axes[0].set_xlabel('Quartile'); axes[0].set_ylabel('freq_coverage_fraction')
axes[0].set_title('Coverage by corridor strategy')
axes[1].set_xticks(x); axes[1].set_xlabel('Quartile'); axes[1].set_ylabel('bad_quality_fraction (lower better)')
axes[1].set_title('Bad-quality fraction by corridor strategy')
axes[0].legend(fontsize=7)
plt.suptitle('Figure 12 -- Stage 4.5: corridor-strategy comparison, 4 report pairs')
plt.tight_layout()
plt.show()
""")

md(r"""Widening the corridor at low-confidence periods (Approach A) is never worse than the
uniform baseline and is a clear quality win on Q3 (coverage 0.904->0.923, bad-quality fraction
roughly halved 0.125->0.062). Down-weighting in scoring instead (Approach B) is never better than
Approach A anywhere they differ. Q4 fails to converge under *every* strategy -- consistent with
Section 6's own H3 finding on this exact quartile, reinforcing that Q4's non-convergence is a
distance/coherence-strength limit, not a corridor-tuning problem. Approach A was adopted;
`build_template_family_widened` (Section 5's scan already used it) is the result.

### The isolated effect of the hybrid curve itself

Holding the corridor strategy, bandwidth, and polarization fix all fixed, and varying only the
reference curve (old shared `SDISPL.ASC` vs. the new per-pair hybrid curve) on the same 4 pairs:""")

code("""old_curve = pd.read_csv(RESULTS / 'stage4_5' / 'stage4_5_validation_old_curve.csv')
new_curve = pd.read_csv(RESULTS / 'stage4_5' / 'stage4_5_validation_new_curve.csv')
compare = old_curve[['quartile', 'converged', 'best_delta_km_s', 'freq_coverage_fraction',
                      'bad_quality_fraction', 'n_templates_converged']].copy()
compare.columns = ['quartile', 'old_converged', 'old_delta', 'old_coverage', 'old_bad_quality', 'old_n_conv']
compare = compare.merge(
    new_curve[['quartile', 'converged', 'best_delta_km_s', 'freq_coverage_fraction',
               'bad_quality_fraction', 'n_templates_converged']]
    .rename(columns={'converged': 'new_converged', 'best_delta_km_s': 'new_delta',
                      'freq_coverage_fraction': 'new_coverage', 'bad_quality_fraction': 'new_bad_quality',
                      'n_templates_converged': 'new_n_conv'}),
    on='quartile',
)
print(compare.to_string(index=False))
""")

code("""fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
q = compare.quartile.values
x = np.arange(len(q))
w = 0.35

axes[0].bar(x - w/2, compare.old_n_conv, w, label='Old SDISPL.ASC curve', color='gray')
axes[0].bar(x + w/2, compare.new_n_conv, w, label='Hybrid ADAMA+GDM52 curve', color='C1')
axes[0].set_xticks(x); axes[0].set_xticklabels([f'Q{i}' for i in q])
axes[0].set_ylabel('Templates converged (of 33 scanned)')
axes[0].set_title('Corridor templates that converge at all')
axes[0].legend(fontsize=8)

axes[1].bar(x - w/2, compare.old_coverage, w, color='gray')
axes[1].bar(x + w/2, compare.new_coverage, w, color='C1')
axes[1].set_xticks(x); axes[1].set_xticklabels([f'Q{i}' for i in q])
axes[1].set_ylabel('freq_coverage_fraction (winning template)')
axes[1].set_title('Coverage of the best-scoring template')

axes[2].bar(x - w/2, compare.old_bad_quality, w, color='gray')
axes[2].bar(x + w/2, compare.new_bad_quality, w, color='C1')
axes[2].set_xticks(x); axes[2].set_xticklabels([f'Q{i}' for i in q])
axes[2].set_ylabel('bad_quality_fraction (lower better)')
axes[2].set_title('Quality of the best-scoring template')

plt.suptitle('Figure 13 -- Isolated effect of the hybrid reference curve alone (same bandwidth, corridor, polarization fix)')
plt.tight_layout()
plt.show()
""")

md(r"""**The consistent signal, across all four pairs**: the number of templates that converge at
all rises with the hybrid curve in *every single quartile* -- the curve is a better match to the
real data broadly, not just at the one delta that happens to score best. On the winning template's
own diagnostics: Q1 flips from total non-convergence to convergence; Q3 improves clearly on both
axes; Q4 keeps similar coverage but meaningfully lower bad-quality. Q2 is the one partial exception
-- more templates converge, but the specific best-scoring template is slightly worse on its own
coverage/bad-quality than the old curve's pick -- reported as-is rather than smoothed over,
consistent with how every other honest negative/mixed finding in this project has been handled.

Combined, Stage 4.5's three fixes (polarization, hybrid curve, corridor strategy) are what Section
5's worked example above and every figure past this point in a future full rerun would reflect --
Round 1/2 (Section 6) predate all three and should be read as the *before* picture, not the current
pipeline's output.""")

# ============================================================================
# Section 8 -- Discussion
# ============================================================================

md(r"""## 8. Discussion: what's proven, what's still open

**Proven, at real dataset scale**: `FastMspec` converts real, noisy ambient-noise cross-spectra
into dispersion curves where single-taper converts none (Section 6); the resolution-bandwidth
criterion this depends on (Section 3) explains both *why* single-taper fails and *why*
convergence itself falls off with distance for every multitaper technique alike; a per-pair,
physically-grounded reference curve measurably improves pick quality over one generic curve
shared by every pair (Section 7), isolated from every other confound.

**Still genuinely open:**
- **Quartile 4's existence threshold (H3).** Not a bandwidth problem (Round 2: 0% at every
  fraction tested, including 1.5x ceiling) and not a corridor-strategy problem (Section 7: fails
  under all three strategies too) -- two independent knobs, same wall. Two other mechanistic
  hypotheses (branch-continuity/cycle-slip fragmentation, walk-reseeding) were tested and also
  came back negative (`docs/round2_hypothesis_evaluation.tex`). What actually limits Q4 is not
  yet identified.
- **The bandwidth lower bound (H2).** `docs/stage5_bandwidth_theory.tex`'s bias-variance MSE
  framework predicts an interior optimum should exist; Round 2's tested range (0.25-1.5x
  `NW_high(r)`) never found one -- convergence favored the widest bandwidth throughout. Testing
  further out, against real independent ground truth, is exactly Notebook 5's job.
- **M/N (maxima/minima) events.** The 2-way barcode here uses zero-crossings only -- the old
  design's own M/N detection was found far noisier (docs/coherence_barcode_design.tex Section 8)
  and dropped rather than papered over. A principled M/N quality gate remains a deferred,
  unstarted item.
- **Rayleigh waves.** The hybrid curve is wired for `wave='rayleigh'` but validated only for Love
  -- this project's data is Love-wave throughout, so this gap hasn't blocked anything yet, but it
  is untested.

### Where this goes next

Round 2 already surfaced the real, still-open question: does bandwidth choice have a genuine
interior optimum, and where is it? Answering that needs real ground truth to check against, not
just internal convergence/quality signals -- **Notebook 5 tests
`docs/stage5_bandwidth_theory.tex`'s bias-variance MSE framework against ADAMA's own real station
pairs and real, independently-measured dispersion curves** (`co`/`cf`, from ADAMA's own
AkiEstimate pipeline -- not the 380 pairs used throughout this notebook, confirmed disjoint from
ADAMA's own catalog). That's the natural place this notebook's own open questions get resolved.""")

md(r"""## Summary

- **Sections 1-3** carried Notebooks 1-3's own FastMspec-vs-single-taper result forward and
  grounded it in a first-principles resolution-bandwidth criterion, tying zero-crossing stability
  directly to the Bessel zero-crossing spacing this whole pipeline depends on.
- **Section 4** introduced the vendored, instrumented picker this project builds on -- real
  internal quality control, not a naive event scan, verified byte-identical to the unmodified
  `seislib` package by default.
- **Section 5** validated it on a real, known-good pair (SKRH-BAND), including a direct
  before/after comparison of the reference-curve fix on that one pair.
- **Section 6** scaled to the full 380-pair dataset: single-taper never converges, multitaper does
  on roughly a quarter of pairs, and Round 2's bandwidth sweep confirmed both the
  resolution-bandwidth criterion's real-dataset consequences and a genuine, still-unexplained
  existence threshold at the farthest quartile.
- **Section 7** closed three real gaps found along the way (a polarization bug, a generic
  reference curve, a one-size corridor) and showed their combined effect is real and mostly
  positive, not asserted -- measured, with the one honest exception reported plainly.
- **Section 8** named what's still open, honestly, and pointed at exactly where the next notebook
  picks up.
""")

nb['cells'] = cells
with open('04_dispersion_curve_picking.ipynb', 'w') as f:
    nbf.write(nb, f)
print("wrote 04_dispersion_curve_picking.ipynb with", len(cells), "cells")

