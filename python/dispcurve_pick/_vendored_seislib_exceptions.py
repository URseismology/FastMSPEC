"""Vendored from seislib.exceptions._exceptions (Fabrizio Magrini, seislib==1.2.1) -- exactly the
3 exception classes _vendored_seislib_an_processing.py actually raises/imports. Confirmed
functionally identical to the installed package via a direct source diff at vendoring time
(DispersionCurveException is literally byte-identical; the other two differ only in cosmetic
operator spacing, no logic changed). See NOTES.md
"Why the seislib package dependency was removed" for why this exists as a standalone file instead
of `from seislib.exceptions import ...`: the full seislib package pulls in an unrelated, broken
Cython extension (seislib.tomography._ray_theory) that fails to build on some HPC toolchains
(bluehive's login-node gcc, specifically) -- these 3 small, stable classes don't need any of that.
"""


class DispersionCurveException(Exception):
    """
    Exception raised when a dispersion curve could not be extracted from the
    data
    """

    def __init__(self):
        self.message = 'It was not possible to retrieve a dispersion curve'
        super().__init__(self.message)

    def __str__(self):
        return self.message


class TimeSpanException(Exception):
    """
    Exception raised when no common time span is found in two obspy traces or
    streams.
    """

    def __init__(self, *args, message=None):
        if message is not None:
            self.message = message
        else:
            self.message = 'No common time span found.'
        for arg in args:
            self.message += '\n%s' % arg
        super().__init__(self.message)

    def __str__(self):
        return self.message


class NonFiniteDataException(Exception):
    """
    Exception raised when the data should be strictly finite but contain
    infinite or nan values.
    """

    def __init__(self, *args):
        self.message = 'The data should be strictly finite, but contain either'
        self.message += ' infinite or nan values.'
        for arg in args:
            self.message += '\n%s' % arg
        super().__init__(self.message)

    def __str__(self):
        return self.message
