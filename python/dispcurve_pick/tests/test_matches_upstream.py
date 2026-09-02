"""Translation-fidelity test for the vendored, instrumented picker (Stage 2 verification).

This is the direct analogue of ccf_pipeline's Octave-comparison tests, adapted for a
Python-to-Python "did I change behavior" check rather than a MATLAB-to-Python translation check:
asserts the instrumented `extract_dispcurve` produces byte-identical `(crossings, dispersion_curve)`
output to the unmodified, pip-installed `seislib==1.2.1` package on the same input, both with
`return_diagnostics=False` (the default -- must match trivially) and `return_diagnostics=True`
(crossings/curve must still match exactly; only the extra diagnostics element is new). This test
itself still imports the real `seislib` package (for the comparison), but the picker under test
(`dispcurve_pick.extract_dispcurve`) does not -- it uses this project's own vendored copies of the
handful of seislib utility functions/exceptions it needs, precisely so it doesn't depend on the
full `seislib` package (whose unrelated `tomography` extension fails to build on some HPC
toolchains) -- see `_vendored_seislib_exceptions.py`'s docstring and NOTES.md.

A real-data check against Sayan Swar's own SKRH-BAND result happens separately in Stage 3, once
that data is pulled locally -- this test is deliberately self-contained (synthetic input only) so
it can run without any external data dependency.
"""
import numpy as np
import pytest
from scipy.special import j0

import seislib.an as upstream_seislib
import seislib.exceptions as seislib_exceptions

from dispcurve_pick import _vendored_seislib_exceptions as vendored_exceptions

from dispcurve_pick import extract_dispcurve, PickDiagnostics, DispersionCurveExceptionWithDiagnostics


def _bessel_coherence(freqs, dist_km, c_of_f):
    """Synthetic real-part-of-coherence spectrum from the Aki/Ekstrom Bessel model,
    Re{coherence(f)} = J0(2*pi*f*r/c(f)) -- the same construction used elsewhere in this
    project (nb5_helpers.template_barcode, Notebook 2). No noise: a clean synthetic case is
    exactly what a translation-fidelity check should use (a noisy real case belongs in Stage 3's
    validation against Sayan's own known-good result, not here).
    """
    x = 2 * np.pi * freqs * dist_km / c_of_f(freqs)
    return j0(x)


@pytest.fixture
def synthetic_input():
    dist_km = 300.0
    freqs = np.linspace(0.01, 0.4, 2000)
    c_of_f = lambda f: 3.0 - 1.5 * f  # mild normal dispersion, plausible Love-wave range
    coherence = _bessel_coherence(freqs, dist_km, c_of_f)
    ref_curve = np.column_stack([freqs[::50], c_of_f(freqs[::50])])
    kwargs = dict(
        frequencies=freqs,
        corr_spectrum=coherence,
        interstation_distance=dist_km,
        ref_curve=ref_curve,
        freqmin=0.02,
        freqmax=0.35,
        cmin=1.0,
        cmax=5.0,
        pick_threshold=0,
    )
    return kwargs


def test_default_matches_upstream_exactly(synthetic_input):
    """return_diagnostics=False (the default): must be byte-identical to unmodified upstream."""
    ours = extract_dispcurve(**synthetic_input)
    theirs = upstream_seislib.extract_dispcurve(**synthetic_input)

    assert len(ours) == 2 == len(theirs)
    np.testing.assert_array_equal(ours[0], theirs[0])  # crossings
    np.testing.assert_array_equal(ours[1], theirs[1])  # dispersion_curve


def test_diagnostics_true_curve_still_matches_upstream(synthetic_input):
    """return_diagnostics=True: crossings/curve must still match upstream exactly; only the
    extra diagnostics element is new."""
    ours = extract_dispcurve(**synthetic_input, return_diagnostics=True)
    theirs = upstream_seislib.extract_dispcurve(**synthetic_input)

    assert len(ours) == 3
    np.testing.assert_array_equal(ours[0], theirs[0])
    np.testing.assert_array_equal(ours[1], theirs[1])

    diag = ours[2]
    assert isinstance(diag, PickDiagnostics)
    assert diag.converged is True
    assert 0.0 <= diag.bad_quality_fraction <= 1.0
    assert diag.n_candidate_crossings > 0
    assert diag.n_accepted_picks > 0
    assert 0.0 <= diag.freq_coverage_fraction <= 1.0
    assert diag.mean_amp_ratio > 0


def test_diagnostics_on_failure_carries_diagnostics(synthetic_input):
    """A candidate too weak/short to converge should still raise
    DispersionCurveExceptionWithDiagnostics (a DispersionCurveException subclass) carrying
    whatever diagnostics were captured before the failure, when return_diagnostics=True; and this
    project's own vendored (not upstream's) DispersionCurveException otherwise -- our
    extract_dispcurve no longer imports the real seislib package at all (see
    _vendored_seislib_exceptions.py's docstring: it pulls in an unrelated, broken Cython extension
    on some HPC toolchains), so the class raised is a distinct-but-functionally-identical vendored
    copy, not literally `seislib_exceptions.DispersionCurveException` -- confirmed equal by
    message text, the one thing that actually matters behaviorally, not by `isinstance`."""
    kwargs = dict(synthetic_input)
    # Restrict to a band far too narrow to pass the coverage acceptance test.
    kwargs["freqmin"], kwargs["freqmax"] = 0.20, 0.21

    with pytest.raises(seislib_exceptions.DispersionCurveException) as upstream_excinfo:
        upstream_seislib.extract_dispcurve(**kwargs)

    with pytest.raises(vendored_exceptions.DispersionCurveException) as ours_excinfo:
        extract_dispcurve(**kwargs)  # return_diagnostics=False
    assert str(ours_excinfo.value) == str(upstream_excinfo.value)

    with pytest.raises(DispersionCurveExceptionWithDiagnostics) as excinfo:
        extract_dispcurve(**kwargs, return_diagnostics=True)
    diag = excinfo.value.diagnostics
    assert isinstance(diag, PickDiagnostics)
    assert diag.converged is False
