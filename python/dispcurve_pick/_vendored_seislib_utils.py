"""Vendored from seislib.utils (Fabrizio Magrini, seislib==1.2.1) -- exactly the 3 functions
_vendored_seislib_an_processing.py actually imports, plus the `resample` helper
`adapt_sampling_rate` itself calls. Functionally identical to the installed package -- confirmed
by a direct source diff at vendoring time, which found only cosmetic differences (operator
spacing, some docstring "Notes" sections trimmed); no logic changed. See
_vendored_seislib_exceptions.py's docstring / NOTES.md for why this exists as a standalone file
instead of `from seislib.utils import ...`.
"""
import numpy as np
from obspy import Stream, Trace

from ._vendored_seislib_exceptions import TimeSpanException


def resample(x, fs):
    """Vendored from seislib.utils._utils.resample (not seislib.utils.adapt_sampling_rate's own
    top-level source -- that function just calls this one, defined separately in the same
    upstream module). Lowpass-filters below the target Nyquist frequency, then interpolates."""
    nyquist_f = fs / 2 - (fs / 2) * 0.01
    try:
        x.filter('lowpass', freq=nyquist_f, corners=4, zerophase=True)
    except ValueError:
        pass  # when fs > sampling_rate(x), filtering is not needed
    x.interpolate(sampling_rate=fs, method="weighted_average_slopes")
    return x


def adapt_timespan(st1, st2):
    """
    Slices all traces from the two input streams to the overlapping timerange.
    Then returns a copy of the sliced streams.

    Parameters
    ----------
    st1, st2 : obspy.Stream, obspy.Trace

    Returns
    -------
    st1, st2 : obspy.Stream, obspy.Trace
        Obspy stream or trace depending on the input. The original input is
        not permanently modified (a copy is returned)

    Raises
    ------
    TimeSpanException
        If no overlap is available
    """
    is_trace = False
    if isinstance(st1, Trace) or isinstance(st2, Trace):
        is_trace = True
    st1 = Stream(st1) if isinstance(st1, Trace) else st1
    st2 = Stream(st2) if isinstance(st2, Trace) else st2

    # this has to be done twice, because otherwise there sometimes occurs a 1s timeshift
    for adapt in range(2):
        starttime = max([tr.stats.starttime for tr in st1]
                         + [tr.stats.starttime for tr in st2])
        endtime = min([tr.stats.endtime for tr in st1]
                       + [tr.stats.endtime for tr in st2])
        if starttime >= endtime:
            raise TimeSpanException(st1, st2)

        st1 = st1.slice(starttime, endtime)
        st2 = st2.slice(starttime, endtime)
        for tr in st1:
            tr.stats.starttime = starttime
        for tr in st2:
            tr.stats.starttime = starttime

    return (st1, st2) if not is_trace else (st1[0], st2[0])


def adapt_sampling_rate(st1, st2):
    """
    If the input streams (or traces) have different sampling rates, the one
    characterized by the largest sampling rate is downsampled to the sampling
    rate of the other stream (or trace).

    Parameters
    ----------
    st1, st2 : obspy.Stream, obspy.Trace

    Returns
    -------
    st1, st2 : obspy.Stream, obspy.Trace
        Obspy stream or trace depending on the input. The input is permanently
        modified
    """
    is_trace = False
    if isinstance(st1, Trace) or isinstance(st2, Trace):
        is_trace = True
    st1 = Stream(st1) if isinstance(st1, Trace) else st1
    st2 = Stream(st2) if isinstance(st2, Trace) else st2
    fs1, fs2 = st1[0].stats.sampling_rate, st2[0].stats.sampling_rate
    if fs1 < fs2:
        st2 = resample(st2, fs1)
    elif fs2 < fs1:
        st1 = resample(st1, fs2)
    return (st1, st2) if not is_trace else (st1[0], st2[0])


def running_mean(x, N):
    """ Moving average

    Parameters
    ----------
    x : ndarray of shape (m,)
        Data vector

    N : int
        Controls the extent of the smoothing (larger values correspond to larger
        smoothing)

    Returns
    -------
    runmean : ndarray of shape (m,)
        Smoothed input
    """
    if N % 2 == 0:
        N += 1
    idx0 = int((N - 1) / 2)
    runmean = np.zeros(len(x))
    cumsum = np.cumsum(np.insert(x, 0, 0))
    runmean[idx0:-idx0] = (cumsum[N:] - cumsum[:-N]) / N
    for i in range(idx0):
        runmean[i] = np.mean(x[:2 * i + 1])
        runmean[-i - 1] = np.mean(x[-2 * i - 1:])
    return runmean
