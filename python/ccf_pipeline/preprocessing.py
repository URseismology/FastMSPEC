"""Translation of ccf_preprocess_filter_data.m's time-domain stages
(detrend, taper, prefilter) plus their _3dim helpers
(ccf_detrend_3dim.m, ccf_cos_taper_3dim.m, ccf_butterfilt_3dim.m).

Scope note: ccf_preprocess_filter_data.m also has IsMultiTaper/IsFTN/
IsOBN/IsSpecWhiten stages that produce frequency-domain output for a
*different* downstream pathway (the outer ccf_compute_crosscorr_Z.m
dispatcher's first branch, which treats data as already-FFT'd). The real
production config found in this codebase (a2_ccf_run_crosscorr_T_mdg.m)
uses IsMspec=1 with all four of those flags off, i.e. only detrend/taper/
prefilter actually run before handing off time-domain data to the
ccf_compute_crosscorr_mtc_* functions this package already translates.
Those four frequency-domain-producing stages are therefore out of scope
here (not needed for the IsMspec pathway); see NOTES.md.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import detrend as _scipy_detrend
from scipy.signal import butter, filtfilt


def ccf_detrend_3dim(data: np.ndarray) -> np.ndarray:
    """Linear detrend of each (day, window) trace along the sample axis.
    MATLAB's detrend() operates column-wise on the (samples, traces)
    reshape; scipy.signal.detrend along axis=2 of the original 3-D array
    is exactly equivalent (detrending is independent per-trace regardless
    of how traces are ordered/reshaped).
    """
    return _scipy_detrend(data, axis=2, type="linear")


def ccf_cos_taper_3dim(data: np.ndarray) -> np.ndarray:
    """5%-of-length half-cosine taper on each end, shared across all
    (day, window) traces, applied along the sample axis (dim 2).
    """
    datalen = data.shape[2]
    m = int(np.floor((datalen * 5 / 100) / 2 + 0.5))
    j = np.arange(1, datalen + 1)  # 1-based, matching the .m source exactly

    tapered = np.zeros(datalen)
    # Sequential if/elseif/elseif in the .m source, translated as exclusive masks
    # in the same precedence order (not independently-derived conditions).
    rising = j <= m + 1
    flat = (~rising) & (j < datalen - m - 1)
    falling = (~rising) & (~flat)

    tapered[rising] = 0.5 * (1 - np.cos(j[rising] * np.pi / (m + 1)))
    tapered[flat] = 1.0
    tapered[falling] = 0.5 * (1 - np.cos((datalen - j[falling]) * np.pi / (m + 1)))

    return data * tapered.reshape(1, 1, datalen)


def ccf_butterfilt_3dim(data: np.ndarray, frange_prefilt, dt: float) -> np.ndarray:
    """2nd-order Butterworth bandpass + zero-phase (filtfilt) filtering
    along the sample axis. MATLAB's `butter(2, frange*2*dt)` maps directly
    onto scipy.signal.butter with the same Wn-normalized-to-Nyquist
    convention -- confirmed via freqz self-check (-3dB points land within
    ~0.001 of the requested cutoffs, passband gain ~1.0, and impulse
    response is exactly symmetric i.e. zero-phase).

    CAVEAT (see NOTES.md): the original .m source's FiltFiltM.m is a
    custom, "rewritten from scratch" implementation of Gustafsson's (1996)
    initial-condition method, not the textbook version. Neither scipy's
    default (reflect-padding) nor its `method='gust'` filtfilt variant
    matches it closely even away from the signal edges (checked against a
    real Octave run: ~1-2% relative error in the middle of a 500-sample
    test signal, worse near the edges) -- this looks like a genuine
    algorithmic difference in FiltFiltM's specific matrix formulation, not
    just an edge-padding convention mismatch. Using scipy's standard
    filtfilt (default method) here as the best available idiomatic Python
    equivalent. Low priority to chase further: every real production
    config found in this codebase has `IsPrefilter=0`, so this function
    is not on the exercised path.
    """
    wn = np.asarray(frange_prefilt) * 2 * dt
    b, a = butter(2, wn, btype="bandpass")
    return filtfilt(b, a, data, axis=2)
