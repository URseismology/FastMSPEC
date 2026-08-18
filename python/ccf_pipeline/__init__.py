"""Python port of the CCF multitaper cross-correlation pipeline
(lib/ccf_compute_crosscorr_mtc_Z.m / _T.m), built on top of the already
translated thomson_multitaper.FastMultitaper library.

See docs/plan_ccf_mtc_translation.md and NOTES.md in this package for the
phased approach and translation caveats.
"""
from .fast_cross_spectrum import fast_spectrum_batch
from .classical_cross_spectrum import classical_spectrum_batch
from .crosscorr_mtc import (
    compute_crosscorr_mtc_fastmspec,
    compute_crosscorr_mtc_mspec,
    compute_crosscorr_mtc_mspecbestk,
    CrosscorrResult,
)
from .dispatch import compute_crosscorr, FilterConfig
from .preprocessing import ccf_detrend_3dim, ccf_cos_taper_3dim, ccf_butterfilt_3dim

__all__ = [
    "fast_spectrum_batch",
    "classical_spectrum_batch",
    "compute_crosscorr_mtc_fastmspec",
    "compute_crosscorr_mtc_mspec",
    "compute_crosscorr_mtc_mspecbestk",
    "CrosscorrResult",
    "compute_crosscorr",
    "FilterConfig",
    "ccf_detrend_3dim",
    "ccf_cos_taper_3dim",
    "ccf_butterfilt_3dim",
]
