#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve St-Venant's problem in 1D by calling a Riemann solver written in C.
See:
https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.python-as-glue.html#index-3
"""

from ctypes import cdll, c_double
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from timeit import default_timer
from riemann import riemann_python
import functools

hL, uL = 2., 0.
hR, uR = 1., 0.

wL = np.array([hL, hL*uL])
wR = np.array([hR, hR*uR])

nx = 1001

xmin = -9.
xmax = 9.
t = 1.0
xi = np.linspace(xmin, xmax, num=nx)/t


def timer(func):
    """Decorator that display the function execution time"""
    @functools.wraps(func)
    def func_call(*args, **kwargs):
        start_time = default_timer()  # Start timer
        result = func(*args, **kwargs)
        elapsed_time = default_timer() - start_time
        print("Time to run {} with {} [s]: {:f}".format(func.__name__,
                                                        args[0].__name__,
                                                        elapsed_time))
        return result
    return func_call


# Load the C-library compiled with "make"
libc = cdll.LoadLibrary("./libstvenant_c.so")
# Explicit the C-function argument types
libc.riemann.argtypes = [np.ctypeslib.ndpointer(float, ndim=1,  # wL
                                                flags='aligned'),
                         np.ctypeslib.ndpointer(float, ndim=1,  # wR
                                                flags='aligned'),
                         c_double,  # xi_j
                         np.ctypeslib.ndpointer(float, ndim=1,  # w
                                                flags='aligned, contiguous,'
                                                      'writeable')]


def riemann_C(wL, wR, xi_j):
    """A wrapper to the C-function libc.riemann"""

    wL = np.require(wL, float, ['ALIGNED'])
    wR = np.require(wR, float, ['ALIGNED'])
    wi = np.zeros_like(wL)
    libc.riemann(wL, wR, xi_j, wi)
    return wi


@timer
def riemann_loop(riemann_func, wL, wR, xi):
    """Loop over xi to return w as a numpy array of size nx using the
    Riemann solver provided by the riemann_func function"""
    return np.array([riemann_func(wL, wR, xi_j) for xi_j in xi])


def load_file(filename="plotriem"):
    """Reads file and return a (x, h, u) tuple"""
    data = pd.read_csv(filename, delim_whitespace=True, header=None).values
    return data[:, 0], data[:, 1], data[:, 2]


def stvenant(plot_file=False):
    """main function that loops over x and plots the results"""

    w_py = riemann_loop(riemann_python, wL, wR, xi)
    w_C = riemann_loop(riemann_C, wL, wR, xi)

    plt.plot(xi, w_py[:, 0], label="h_py")
    plt.plot(xi, w_py[:, 1]/w_py[:, 0], label="u_py")
    plt.plot(xi, w_C[:, 0], label="h_C")
    plt.plot(xi, w_C[:, 1]/w_C[:, 0], label="u_C")
    plt.xlabel(r'$\xi$')

    if plot_file:
        filename = "plotriem"
        xi_f, h_f, u_f = load_file(filename)
        plt.plot(xi_f, h_f, 'g+', label="h_"+filename)
        plt.plot(xi_f, u_f, 'k+', label="u_"+filename)

    plt.legend()
    plt.show()


if __name__ == '__main__':
    stvenant(plot_file=True)
