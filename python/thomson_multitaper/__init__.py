"""Python translation of ThomsonsMethodRevisitedExperiments/ (Santhosh Karnik's
Fast Slepian Transform toolbox, as vendored under
PRJ_SPAC/codes/test/matlab/functions/).

See NOTES.md for translation caveats and unresolved items.
"""
from .fast_multitaper import FastMultitaper
from .first_n_lambda_dpss import first_n_lambda_dpss
from .multitaper import Multitaper
from .multitaper_adaptive import MultitaperAdaptive
from ._dpss_transition import transition_dpss, transition_dpss_modif
from ._tridiagonal import tridieig, tridisolve
from ._utils import datawrap, dpss

__all__ = [
    "FastMultitaper",
    "Multitaper",
    "MultitaperAdaptive",
    "first_n_lambda_dpss",
    "transition_dpss",
    "transition_dpss_modif",
    "tridieig",
    "tridisolve",
    "datawrap",
    "dpss",
]
