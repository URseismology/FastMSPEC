"""Diagnostics types for the instrumented seislib picker.

See NOTES.md for the full provenance/instrumentation story. In short: seislib's own
`extract_dispcurve` (_vendored_seislib_an_processing.py) computes a rich set of internal
quality signals while picking a dispersion curve -- a per-crossing `bad_quality` flag, the
amplitude ratio behind each accepted pick, and the frequency-coverage fraction its own
acceptance test already checks -- but discards all of it, returning only the final curve or
raising an exception. `PickDiagnostics` is what the instrumented picker surfaces instead, when
called with `return_diagnostics=True`.
"""
from __future__ import annotations

from dataclasses import dataclass

from ._vendored_seislib_exceptions import DispersionCurveException


@dataclass
class PickDiagnostics:
    """Quality signals seislib's own picker already computes internally, surfaced explicitly.

    converged: whether picking succeeded (a curve was returned) or the coverage acceptance
        test rejected it (DispersionCurveExceptionWithDiagnostics was raised instead).
    bad_quality_fraction: fraction of candidate zero crossings seislib's own 3-criterion gate
        (peak-ratio, envelope-relative-amplitude, spacing-vs-reference) flagged as bad quality,
        before picking even started. Lower is better.
    n_candidate_crossings: total zero crossings considered (bad-quality and good, combined).
    n_accepted_picks: picks that survived the kernel-density picking loop's own local
        amplitude-ratio/cycle-jump gates, before final smoothing/coverage filtering.
    freq_coverage_fraction: fraction of the log-spaced frequency grid the final smoothed curve
        actually spans -- the same quantity seislib's own acceptance test compares against 1/5.
    mean_amp_ratio: mean, across accepted picks, of the local kernel-density
        maximum-to-flanking-minimum amplitude ratio (`maxamp/minamp`) that gated each pick's
        acceptance (`pick_threshold` in the original function). Higher means each pick sat on a
        more sharply-defined, less ambiguous kernel-density peak. NaN if no picks were accepted.
    """

    converged: bool
    bad_quality_fraction: float
    n_candidate_crossings: int
    n_accepted_picks: int
    freq_coverage_fraction: float
    mean_amp_ratio: float


class DispersionCurveExceptionWithDiagnostics(DispersionCurveException):
    """Raised instead of the plain DispersionCurveException when return_diagnostics=True and
    picking fails -- carries a `.diagnostics` attribute so a failed pick's quality signals
    (which crossings existed, how bad they were) aren't lost along with the exception. This is
    the more common case for a large sweep: many pairs may not converge, and knowing *why* (a
    high bad-quality fraction vs. simply too few candidate crossings) is itself useful signal.
    """

    def __init__(self, diagnostics: PickDiagnostics):
        self.diagnostics = diagnostics
        # DispersionCurveException.__init__ takes no arguments (a fixed message) -- call it as
        # upstream defines it, then override self.message with a richer, diagnostics-derived one
        # (its __str__ just returns self.message, so this is a clean override, not a workaround).
        super().__init__()
        self.message = (
            f"Dispersion curve picking did not converge (bad_quality_fraction="
            f"{diagnostics.bad_quality_fraction:.2f}, n_candidate_crossings="
            f"{diagnostics.n_candidate_crossings}, freq_coverage_fraction="
            f"{diagnostics.freq_coverage_fraction:.2f})"
        )
