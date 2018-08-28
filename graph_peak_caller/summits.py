import numpy as np


# Modified from http://www.scipy.org/Cookbook/SavitzkyGolay
# positive window_size not enforced anymore
# needs sane input paramters, window size > 4
# switched to double precision for internal accuracy

def find_summits(signal, window_size=51):
    """return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified
    """
    window_size = (window_size//2)*2+1  # to make a even number
    sg = savitzky_golay_order2(signal, window_size, deriv=1)
    return np.where(np.diff(np.sign(sg)) <= -1)[0]


def savitzky_golay_order2(signal, window_size, deriv=0):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    if window_size % 2 != 1:
        window_size += 1
    half_window = (window_size - 1) // 2

    # precompute coefficients
    b = np.mat([[1, k, k**2] for k in range(-half_window, half_window+1)],
               dtype='int64')
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
    lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
    signal = np.concatenate((firstvals, signal, lastvals))
    ret = np.convolve(m[::-1], signal.astype('float64'), mode='valid').astype('float32')
    return ret
