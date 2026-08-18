"""Validation for ccf_pipeline.prepare_data (Phase 5).

No real SAC files or Octave fixtures available while writing this (see
NOTES.md) -- verification here is self-consistency against the .m
source's own windowing formula (nwin, win_length) and obspy-synthesized
Trace objects, not a real Octave/data cross-check. Re-verify against real
data and/or Octave once SAC files are available.
"""
import numpy as np
from obspy import Trace, UTCDateTime

from ccf_pipeline.prepare_data import build_windows, validate_pair, station_distance_km


def _make_trace(data, dt=1.0, stla=10.0, stlo=20.0, stel=0.0, t0=None):
    t0 = t0 or UTCDateTime(2024, 1, 1)
    return Trace(
        data=np.asarray(data, dtype=float),
        header={"delta": dt, "starttime": t0, "npts": len(data), "sac": {"stla": stla, "stlo": stlo, "stel": stel}},
    )


def test_validate_pair_accepts_good_data():
    npts = 86500
    tr1 = _make_trace(np.random.randn(npts))
    tr2 = _make_trace(np.random.randn(npts))
    assert validate_pair(tr1, tr2, dt=1.0) is None


def test_validate_pair_rejects_all_zero():
    npts = 86500
    tr1 = _make_trace(np.zeros(npts))
    tr2 = _make_trace(np.random.randn(npts))
    assert validate_pair(tr1, tr2, dt=1.0) == "all zeros"


def test_validate_pair_rejects_too_short():
    tr1 = _make_trace(np.random.randn(100))
    tr2 = _make_trace(np.random.randn(100))
    reason = validate_pair(tr1, tr2, dt=1.0)
    assert reason is not None and "too short" in reason


def test_validate_pair_errors_on_delta_mismatch():
    tr1 = _make_trace(np.random.randn(86500), dt=1.0)
    tr2 = _make_trace(np.random.randn(86500), dt=2.0)
    try:
        validate_pair(tr1, tr2, dt=1.0)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_build_windows_shape_matches_matlab_formula():
    npts = 86500
    dt = 1.0
    tr1 = _make_trace(np.random.randn(npts), dt=dt)
    tr2 = _make_trace(np.random.randn(npts), dt=dt)

    winlength_hours = 3
    nstart_sec = 50
    s1w, s2w = build_windows(tr1, tr2, winlength_hours, nstart_sec, dt)

    hour_length = winlength_hours
    nwin_expected = int(np.floor(24 / hour_length)) * 2 - 1
    win_length = int(hour_length * 3600 / dt)
    last_pt = win_length * 0.5 * (nwin_expected - 1) + 1 + (nstart_sec / dt) + win_length
    if last_pt < npts:
        nwin_expected += 1

    assert s1w.shape == (nwin_expected, win_length + 1)
    assert s2w.shape == (nwin_expected, win_length + 1)


def test_build_windows_boundaries_align_with_formula():
    """Uses a ramp signal (data[k]=k) so window contents directly reveal
    which sample indices were cut, independent of interpolation subtleties.
    """
    npts = 86500
    dt = 1.0
    tr1 = _make_trace(np.arange(npts, dtype=float), dt=dt)
    tr2 = _make_trace(np.arange(npts, dtype=float), dt=dt)

    winlength_hours = 3
    nstart_sec = 50
    s1w, _ = build_windows(tr1, tr2, winlength_hours, nstart_sec, dt)

    win_length = int(winlength_hours * 3600 / dt)
    # window 1 (0-indexed 0): pts_begin (1-based) = 1 + nstart_sec/dt
    assert s1w[0, 0] == 1 + nstart_sec / dt
    # window 2 (0-indexed 1): pts_begin = win_length*0.5*1 + 1 + nstart_sec/dt
    assert s1w[1, 0] == win_length * 0.5 + 1 + nstart_sec / dt


def test_station_distance_reasonable():
    # ~0.5 deg lat/lon apart near the equator is roughly 50-80 km
    d = station_distance_km(10.0, 20.0, 10.5, 20.5)
    assert 50 < d < 90


if __name__ == "__main__":
    import sys

    tests = [
        test_validate_pair_accepts_good_data,
        test_validate_pair_rejects_all_zero,
        test_validate_pair_rejects_too_short,
        test_validate_pair_errors_on_delta_mismatch,
        test_build_windows_shape_matches_matlab_formula,
        test_build_windows_boundaries_align_with_formula,
        test_station_distance_reasonable,
    ]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {t.__name__}: {exc}")
    print()
    print("ALL PASS" if not failures else f"{failures} FAILURE(S)")
    sys.exit(1 if failures else 0)
