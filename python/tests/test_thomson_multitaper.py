"""Validation suite for the ThomsonsMethodRevisitedExperiments Python translation.

Cross-checks against independent references (numpy dense eigensolvers,
scipy's own DPSS implementation) rather than against the MATLAB source
directly, since no MATLAB/Octave was available on the translating machine.
Run with: python3 -m pytest tests/ -v  (or plain `python3 tests/test_thomson_multitaper.py`)
"""
import numpy as np

from thomson_multitaper import (
    FastMultitaper,
    Multitaper,
    MultitaperAdaptive,
    dpss,
    first_n_lambda_dpss,
    tridieig,
    tridisolve,
    transition_dpss,
)


def test_tridieig_matches_dense_eigensolver():
    rng = np.random.default_rng(0)
    n = 8
    diag = rng.standard_normal(n)
    off = rng.standard_normal(n - 1)
    a = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    ref = np.sort(np.linalg.eigvalsh(a))
    mine = tridieig(diag, np.concatenate([[0.0], off]), 1, n)
    assert np.allclose(ref, mine, atol=1e-10)


def test_tridisolve_matches_dense_solve():
    rng = np.random.default_rng(1)
    n = 10
    diag = rng.standard_normal(n) + 5
    off = rng.standard_normal(n - 1)
    a = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    b = rng.standard_normal(n)
    ref = np.linalg.solve(a, b)
    mine = tridisolve(off, diag, b)
    assert np.allclose(ref, mine, atol=1e-10)


def test_transition_dpss_matches_scipy():
    for n, w in [(64, 4 / 64), (128, 6 / 128), (200, 3 / 200)]:
        _, lam, _, lo, up = transition_dpss(n, w)
        _, lam_ref = dpss(n, n * w, up)
        assert np.allclose(lam, lam_ref[lo - 1 : up], atol=1e-6)


def test_first_n_lambda_dpss_matches_scipy():
    for n, w, first_n in [(256, 4 / 256, 8), (100, 5 / 100, 15), (64, 2 / 64, 6)]:
        lam = first_n_lambda_dpss(n, w, first_n)
        _, lam_ref = dpss(n, n * w, first_n)
        assert np.allclose(lam, lam_ref, atol=1e-5)


def _synthetic_tone(n, f0, seed=2, noise=0.3):
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    return np.sin(2 * np.pi * f0 * t) + noise * rng.standard_normal(n)


def test_multitaper_recovers_tone_peak():
    n, w, f0 = 256, 4 / 256, 0.08
    x = _synthetic_tone(n, f0)
    freqs = np.fft.fftfreq(n, d=1)
    mt = Multitaper(n, w, cutoff=6)
    y = mt.spectral_estimate(x, m=n)
    peak = abs(freqs[np.argmax(y[: n // 2])])
    assert abs(peak - f0) < 0.01


def test_multitaper_adaptive_recovers_tone_peak():
    n, w, f0 = 256, 4 / 256, 0.08
    x = _synthetic_tone(n, f0)
    freqs = np.fft.fftfreq(n, d=1)
    mta = MultitaperAdaptive(n, w, cutoff=6)
    sf, _, _, _ = mta.spectral_estimate(x, m=n)
    peak = abs(freqs[np.argmax(sf[: n // 2])])
    assert abs(peak - f0) < 0.01


def test_fast_multitaper_matches_exact_multitaper():
    n, w = 256, 4 / 256
    x = _synthetic_tone(n, 0.08)
    mt = Multitaper(n, w, cutoff=6)
    fmt = FastMultitaper(n, w, cutoff=6)
    y = mt.spectral_estimate(x, m=n)
    z = fmt.spectral_estimate(x, m=n)
    assert np.linalg.norm(y - z) / np.linalg.norm(y) < 1e-6


def test_decimated_output_path_m_less_than_n():
    n, w, m = 256, 4 / 256, 64
    x = _synthetic_tone(n, 0.08)
    for estimator in (Multitaper(n, w, cutoff=6), FastMultitaper(n, w, cutoff=6)):
        y = estimator.spectral_estimate(x, m=m)
        assert y.shape == (m,)
        assert np.all(np.isfinite(y))


def test_fractional_cutoff_mode_runs():
    n, w = 256, 4 / 256
    x = _synthetic_tone(n, 0.08)
    for cls in (Multitaper, FastMultitaper):
        est = cls(n, w, cutoff=0.05)
        y = est.spectral_estimate(x, m=n)
        assert np.all(np.isfinite(y))


if __name__ == "__main__":
    import sys

    tests = [obj for name, obj in list(globals().items()) if name.startswith("test_")]
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
