"""Validation for ccf_pipeline.dispatch (Phase 6).

The IsMspec routing is verified against a real Octave run of the actual
ccf_compute_crosscorr_Z.m (the dispatcher itself, not just _mtc_Z.m
directly) -- confirms the wiring, not just the underlying computation
already verified in test_crosscorr_mtc.py. The other two branches
(plain-fft, already-frequency-domain) are smoke-tested only: they are not
the production path (every real config found sets IsMspec=1), and were
not built with an Octave cross-check the way the IsMspec path was.
"""
import numpy as np
import scipy.io as sio
import pytest

from ccf_pipeline import compute_crosscorr, FilterConfig

FIXTURE = "../octave_verify_ccf/synth_medium.mat"
DISPATCH_FIXTURE = "../octave_verify_ccf/dispatch_fastmspec_out.mat"


def _load(path):
    try:
        return sio.loadmat(path)
    except FileNotFoundError:
        pytest.skip(f"Octave fixture {path} not present")


def test_ismspec_routing_matches_octave_dispatcher():
    fix = _load(FIXTURE)
    s1, s2 = fix["S1"], fix["S2"]
    mat = _load(DISPATCH_FIXTURE)
    coh_octave = mat["coh_sum"].squeeze()

    config = FilterConfig(dt=1.0, is_mspec=True, technique="FastMspec", wband=0.05, cutoff=0.95, epsilon=0.05)
    result = compute_crosscorr(s1, s2, config)

    n_samples = s1.shape[2]
    n_onesided = n_samples // 2 + 1
    rel = np.linalg.norm(result.coh_sum[:n_onesided] - coh_octave[:n_onesided]) / np.linalg.norm(
        coh_octave[:n_onesided]
    )
    assert rel < 1e-10
    assert result.taper_size == 4
    assert result.coh_num == 6


def test_plain_fft_branch_runs():
    rng = np.random.default_rng(0)
    s1, s2 = rng.standard_normal((3, 2, 32)), rng.standard_normal((3, 2, 32))
    config = FilterConfig(dt=1.0)  # all flags off -> plain fft branch
    coh_trace, coh_num = compute_crosscorr(s1, s2, config)
    assert coh_trace.shape == s1.shape
    assert coh_num == 6
    assert np.all(np.isfinite(coh_trace))


def test_already_freq_domain_branch_runs():
    rng = np.random.default_rng(0)
    s1, s2 = rng.standard_normal((3, 2, 32)), rng.standard_normal((3, 2, 32))
    config = FilterConfig(dt=1.0, is_multitaper=True)
    coh_trace, coh_num = compute_crosscorr(s1, s2, config)
    assert coh_trace.shape == s1.shape
    assert coh_num == 6
    assert np.all(np.isfinite(coh_trace))


if __name__ == "__main__":
    import sys

    tests = [
        test_ismspec_routing_matches_octave_dispatcher,
        test_plain_fft_branch_runs,
        test_already_freq_domain_branch_runs,
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
