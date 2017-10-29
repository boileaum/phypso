cimport cstvenant
import numpy as np
from cpython cimport array


def riemann_cython(wL, wR, xi_j):
    """A wrapper to the C-function riemann"""

    cdef array.array wL_arr = array.array('d', wL)
    cdef array.array wR_arr = array.array('d', wR)

    wi = np.zeros_like(wL)
    cdef array.array wi_arr = array.array('d', wi)

    cstvenant.riemann(wL_arr.data.as_doubles, wR_arr.data.as_doubles, xi_j,
                      wi_arr.data.as_doubles)
    return wi_arr
