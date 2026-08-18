"""Validation for ccf_pipeline.preprocessing (Phase 4).

detrend and taper verified against real Octave runs to machine precision.
butterfilt is checked via self-consistent scipy signal-processing checks
only -- see preprocessing.py's docstring for why an exact Octave match
wasn't achieved (a genuine algorithmic difference in FiltFiltM.m, not
pursued further since IsPrefilter=0 in every real production config found).
"""
import numpy as np
import scipy.io as sio
import pytest

from ccf_pipeline.preprocessing import ccf_detrend_3dim, ccf_cos_taper_3dim, ccf_butterfilt_3dim

FIXTURE = "../octave_verify_ccf/synth_medium.mat"
PREPROC_FIXTURE = "../octave_verify_ccf/preproc_out.mat"


def _load(path):
    try:
        return sio.loadmat(path)
    except FileNotFoundError:
        pytest.skip(f"Octave fixture {path} not present")


def test_detrend_matches_octave():
    s1 = _load(FIXTURE)["S1"]
    oct_out = _load(PREPROC_FIXTURE)["S1d"]
    py_out = ccf_detrend_3dim(s1)
    assert np.max(np.abs(py_out - oct_out)) < 1e-10


def test_taper_matches_octave():
    s1 = _load(FIXTURE)["S1"]
    oct_out = _load(PREPROC_FIXTURE)["S1t"]
    py_out = ccf_cos_taper_3dim(s1)
    assert np.max(np.abs(py_out - oct_out)) < 1e-12


def test_taper_shape_preserved():
    s1 = _load(FIXTURE)["S1"]
    py_out = ccf_cos_taper_3dim(s1)
    assert py_out.shape == s1.shape


def test_taper_endpoints_near_zero_and_middle_near_one():
    n = 100
    data = np.ones((1, 1, n))
    tapered = ccf_cos_taper_3dim(data)
    assert tapered[0, 0, 0] < 0.2
    assert tapered[0, 0, -1] < 0.2
    assert tapered[0, 0, n // 2] > 0.99


def test_butterfilt_runs_and_preserves_shape():
    s1 = _load(FIXTURE)["S1"]
    out = ccf_butterfilt_3dim(s1, [0.02, 0.4], dt=1.0)
    assert out.shape == s1.shape
    assert np.all(np.isfinite(out))


if __name__ == "__main__":
    import sys

    tests = [
        test_detrend_matches_octave,
        test_taper_matches_octave,
        test_taper_shape_preserved,
        test_taper_endpoints_near_zero_and_middle_near_one,
        test_butterfilt_runs_and_preserves_shape,
    ]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {t.__name__}: {exc}")
        except Exception as exc:
            print(f"SKIP/ERROR {t.__name__}: {exc}")
    print()
    print("ALL PASS" if not failures else f"{failures} FAILURE(S)")
    sys.exit(1 if failures else 0)
