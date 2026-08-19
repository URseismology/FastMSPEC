"""New (not translated from MATLAB) helper code for Notebook 3. Two kinds:

1. `calc_snr_onesided` -- a direct, faithful port of the real
   `calc_SNR_onesided.m` (found in `sayan-swar-translation/.../functions/`)
   used by the actual `ap_ccf_compute_snr.m` script that produced the SNR
   numbers in Sayan's report. This is a genuine translation (kept here
   rather than in `ccf_pipeline/` since it has not been Octave-verified --
   see the module docstring caveats below), reproducing its exact,
   slightly unusual calling convention (negative group-velocity arguments)
   because that convention is what actually produced the report's numbers.

2. Rotation-to-transverse and NLNM-synthetic helpers -- genuinely new code
   for this notebook, not present anywhere in the MATLAB source (no
   N/E -> R/T rotation exists in this codebase at all; see
   python/ccf_pipeline/NOTES.md's "What's genuinely unverified" section).
   Built using obspy's standard, independently-verified rotation function
   (`obspy.signal.rotate.rotate_ne_rt`) rather than a from-scratch
   reimplementation, and obspy's built-in Peterson (1993) NLNM model
   (`obspy.signal.spectral_estimation.get_nlnm`) rather than hand-copying
   the model's published coefficients.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from obspy import read
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
from obspy.signal.spectral_estimation import get_nlnm
from scipy.special import j0


def calc_snr_onesided(ccf: np.ndarray, grv_min: float, grv_max: float, r: float, dt: float):
    """Direct port of calc_SNR_onesided.m (JBR 7/8/18). `ccf` is a one-sided
    or two-sided frequency-domain coherency array (matches whatever ifft(ccf)
    expects in the .m source -- a full-length array, used as-is).

    Note the call convention used by ap_ccf_compute_snr.m itself passes
    grv_min=-2, grv_max=-5 (negative group velocities) -- this looks
    unusual, but the function's own win_min<15/win_max<50 fallback clamps
    make it well-defined (collapsing to a fixed near-zero-lag causal
    window of ~(0, 50) seconds for interstation distances in the range
    used here). Reproduced exactly since matching that behavior is what
    makes this comparable to the report's own reported SNR numbers.
    """
    snrdata = np.real(np.fft.ifft(ccf))
    snrdata = np.fft.fftshift(snrdata)
    nn = len(snrdata)
    lag = (np.arange(nn) - np.floor(nn / 2)) * dt

    win_min = r / grv_max
    win_max = r / grv_min
    if win_min < 15:
        win_min = 0
    if win_max < 50:
        win_max = 50

    signal_ind = (lag > win_min) & (lag < win_max)
    signal_amp = np.sum(snrdata[signal_ind] ** 2) / np.count_nonzero(signal_ind)
    noise_amp = np.sum(snrdata[~signal_ind] ** 2) / np.count_nonzero(~signal_ind)
    snr = 10 * np.log10(signal_amp / noise_amp)
    return snr, signal_ind


def load_orientation(metadata_dir: Path) -> dict[str, float]:
    df = pd.read_csv(Path(metadata_dir) / "orientation.csv", header=None, names=["sta", "az"])
    return {row.sta.split("-")[-1]: row.az for row in df.itertuples()}


def _correct_sensor_orientation(n: np.ndarray, e: np.ndarray, off_deg: float):
    """Rotates a (possibly misoriented) horizontal N/E pair to true N/E,
    given the sensor's N-component azimuth offset from true north
    (orientation.csv's convention)."""
    off = np.deg2rad(off_deg)
    true_n = n * np.cos(off) - e * np.sin(off)
    true_e = n * np.sin(off) + e * np.cos(off)
    return true_n, true_e


def find_ne_day_pairs(datadir: Path, sta1: str, sta2: str, band: str = "LH") -> list[tuple[str, str, str, str, str]]:
    """Like ccf_pipeline.prepare_data.find_day_pairs, but for E/N component
    pairs at both stations (four files per day, all four required). `band`
    is the SEED channel-band-plus-instrument prefix (e.g. "LH" for MTAN/RUNG's
    long-period data, "BH" for SA53/SA58's broadband data) -- component
    naming isn't consistent across the stations in this dataset."""
    dir1 = Path(datadir) / sta1
    comp_e, comp_n = f"{band}E", f"{band}N"
    pairs = []
    for f1e in sorted(dir1.glob(f"{sta1}.*.{comp_e}.sac")):
        suffix = f1e.name[len(sta1):]  # e.g. ".1994.146.00.00.00.LHE.sac"
        day_id = ".".join(f1e.name.split(".")[1:6])
        f1n = dir1 / f"{sta1}{suffix.replace(comp_e, comp_n)}"
        f2e = Path(datadir) / sta2 / f"{sta2}{suffix}"
        f2n = Path(datadir) / sta2 / f"{sta2}{suffix.replace(comp_e, comp_n)}"
        if f1n.exists() and f2e.exists() and f2n.exists():
            pairs.append((str(f1e), str(f1n), str(f2e), str(f2n), day_id))
    return pairs


def prepare_transverse_pair(datadir: Path, metadata_dir: Path, sta1: str, sta2: str,
                             winlength_hours: float, nstart_sec: float, dt: float,
                             min_samples: int = 30000, band: str = "LH"):
    """Builds rotated transverse-component (day, window, sample) arrays for
    a station pair, by: loading E/N at both stations for each day both have
    complete data, correcting sensor orientation (orientation.csv), rotating
    to R/T using the interstation azimuth (obspy's rotate_ne_rt), then
    reusing ccf_pipeline.prepare_data.build_windows unchanged (it only needs
    obspy Trace-like objects with .stats.npts/.delta and .times()/.data) to
    do the actual sliding-window cut-and-resample.

    Returns (t1_data_mat, t2_data_mat, dist_km) or None if no valid days.
    """
    from ccf_pipeline.prepare_data import build_windows, validate_pair, station_distance_km

    orient = load_orientation(metadata_dir)
    pairs = find_ne_day_pairs(datadir, sta1, sta2, band=band)
    if not pairs:
        return None

    t1_days, t2_days = [], []
    dist_km = None
    for f1e, f1n, f2e, f2n, day_id in pairs:
        e1, n1 = read(f1e)[0], read(f1n)[0]
        e2, n2 = read(f2e)[0], read(f2n)[0]
        if validate_pair(e1, e2, dt, min_samples) is not None:
            continue
        if validate_pair(n1, n2, dt, min_samples) is not None:
            continue

        if dist_km is None:
            lat1, lon1 = e1.stats.sac.stla, e1.stats.sac.stlo
            lat2, lon2 = e2.stats.sac.stla, e2.stats.sac.stlo
            dist_km = station_distance_km(lat1, lon1, lat2, lon2)
            _, az12, az21 = gps2dist_azimuth(lat1, lon1, lat2, lon2)

        n1c, e1c = _correct_sensor_orientation(n1.data, e1.data, orient.get(sta1, 0.0))
        n2c, e2c = _correct_sensor_orientation(n2.data, e2.data, orient.get(sta2, 0.0))
        _, t1_data = rotate_ne_rt(n1c, e1c, az12)
        _, t2_data = rotate_ne_rt(n2c, e2c, az21)

        t1 = e1.copy()
        t1.data = t1_data
        t2 = e2.copy()
        t2.data = t2_data

        w1, w2 = build_windows(t1, t2, winlength_hours, nstart_sec, dt)
        t1_days.append(w1)
        t2_days.append(w2)

    if not t1_days:
        return None
    return np.stack(t1_days, axis=0), np.stack(t2_days, axis=0), dist_km


def bessel_fit_quality(freqs: np.ndarray, coh: np.ndarray, dist_km: float,
                        c_grid: np.ndarray, freqmin: float = 0.0, freqmax: float = 1.0):
    """A frequency-domain, model-based alternative to calc_snr_onesided (see
    Notebook 2, Section 1): Aki (1957)'s prediction is Re{coherence(f)} ~
    J0(2*pi*f*r/c) for a diffuse wavefield at a single (or slowly-varying)
    phase velocity c. Grid-searches c to minimize the RMS residual between
    the observed coherence and the Bessel-function prediction, over the
    band [freqmin, freqmax]. A lower best-fit RMS means the observed
    coherence looks more like a physically-plausible diffuse-field spectrum
    (and thus, indirectly, that a phase-velocity fit against it would be
    better-conditioned) -- independent of calc_snr_onesided's fixed-window
    time-domain heuristic, and directly grounded in the Bessel-coherence
    relationship itself rather than an SNR proxy for it.

    Returns (best_c, best_rms, rms_curve) where rms_curve[i] is the RMS
    residual at c_grid[i].
    """
    band = (freqs >= freqmin) & (freqs <= freqmax) & (freqs > 0)
    f_band = freqs[band]
    obs = np.real(coh[band])

    rms_curve = np.empty(len(c_grid))
    for i, c in enumerate(c_grid):
        pred = j0(2 * np.pi * f_band * dist_km / c)
        rms_curve[i] = np.sqrt(np.mean((obs - pred) ** 2))

    best_idx = int(np.argmin(rms_curve))
    return float(c_grid[best_idx]), float(rms_curve[best_idx]), rms_curve


def envelope_conditioned_coherency(coh: np.ndarray, dt: float, threshold_db: float = -12.0):
    """Group-velocity/envelope windowing of a frequency-domain coherency
    spectrum before Bessel-fitting, motivated by (but NOT a literal port
    of -- both papers were inaccessible behind a paywall from this
    environment) standard envelope/group-velocity windowing practice in
    ambient-noise seismology, of the kind used in Hawkins & Sambridge
    (2019, BSSA) and referenced by Xue & Olugboji (2025, JGR-ML/AkiNet).
    See notebooks/README.md's References section for full citations.

    Unlike calc_snr_onesided's FIXED (0, 50) s window (see this notebook's
    earlier callout on why that's a weak heuristic), this computes the
    actual analytic-signal envelope of the time-domain NCF (via Hilbert
    transform) and keeps only the region where the envelope exceeds
    `threshold_db` relative to its own peak -- adapting to each pair's
    actual coherent-arrival structure instead of assuming a fixed lag
    range works for every distance/velocity regime.

    Returns (conditioned_coh, envelope, mask) -- conditioned_coh is the
    re-transformed one-sided-equivalent frequency-domain array (same shape
    as `coh`), for direct use in bessel_fit_quality.
    """
    from scipy.signal import hilbert

    n = len(coh)
    ncf_time = np.real(np.fft.ifft(coh))
    ncf_time = np.fft.fftshift(ncf_time)

    analytic = hilbert(ncf_time)
    envelope = np.abs(analytic)
    peak_db = 20 * np.log10(np.maximum(envelope, 1e-30) / envelope.max())
    mask = peak_db > threshold_db

    conditioned_time = np.fft.ifftshift(ncf_time * mask)
    conditioned_coh = np.fft.fft(conditioned_time)
    return conditioned_coh, envelope, np.fft.ifftshift(mask)


def nlnm_synthetic(n: int, dt: float, seed: int = 0):
    """Generates a length-n synthetic Gaussian random process whose PSD
    follows Peterson (1993)'s New Low Noise Model, using obspy's built-in
    NLNM table (obspy.signal.spectral_estimation.get_nlnm, which returns
    period-in-seconds vs. acceleration-PSD-in-dB, re 1 (m/s^2)^2/Hz).
    Returns (x, freqs, target_psd) where target_psd is on the same
    linear-power scale as a periodogram/multitaper estimate of x.
    """
    periods, db = get_nlnm()
    nlnm_freqs = 1.0 / periods  # ascending period -> descending frequency
    order = np.argsort(nlnm_freqs)
    nlnm_freqs, db = nlnm_freqs[order], db[order]

    freqs = np.fft.rfftfreq(n, d=dt)
    freqs_safe = np.clip(freqs, nlnm_freqs.min(), nlnm_freqs.max())
    db_interp = np.interp(freqs_safe, nlnm_freqs, db)
    target_psd = 10 ** (db_interp / 10)
    target_psd[0] = target_psd[1]  # avoid a zero-frequency singularity

    rng = np.random.default_rng(seed)
    white = rng.standard_normal(n)
    fwhite = np.fft.rfft(white)
    shaped = fwhite * np.sqrt(target_psd * n / (2 * dt))
    x = np.fft.irfft(shaped, n=n)

    full_freqs = np.fft.fftfreq(n, d=dt)
    full_psd = np.interp(np.abs(full_freqs), freqs, target_psd)
    return x, full_freqs, full_psd
