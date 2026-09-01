"""Vendored, instrumented copy of seislib's dispersion-curve picker.

`extract_dispcurve` here is Sayan Swar's own copy of seislib's internal `_an_processing.py`
(itself essentially identical to the pip-installed `seislib==1.2.1` package -- see NOTES.md for
the exact upstream diff), with a small, clearly-marked instrumentation layer added: an optional
`return_diagnostics=True` that surfaces the picker's own already-computed internal quality
signals (bad-quality-crossing fraction, accepted-pick count, frequency-coverage fraction, mean
pick amplitude ratio) as a `PickDiagnostics` object, instead of discarding them. Default behavior
(`return_diagnostics=False`) is byte-identical to calling the unmodified upstream function --
see `tests/test_matches_upstream.py`.

Full design rationale: docs/coherence_barcode_design.tex, Section 8 ("Revision: From Barcode
Matching to Instrumented Reference-Guided Picking").
"""
from ._vendored_seislib_an_processing import extract_dispcurve
from .diagnostics import PickDiagnostics, DispersionCurveExceptionWithDiagnostics

__all__ = [
    "extract_dispcurve",
    "PickDiagnostics",
    "DispersionCurveExceptionWithDiagnostics",
]
