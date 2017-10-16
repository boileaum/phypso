#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve St-Venant's problem in 1D by calling a Riemann solver written in C
(See https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.python-as-glue.html#index-3)
"""

from ctypes import cdll, c_double
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from timeit import default_timer


hL, uL = 2., 0.
hR, uR = 1., 0.

wL = np.array([hL, hL*uL])
wR = np.array([hR, hR*uR])

nx = 1000

w = np.zeros((nx, 2))

xmin = -9.
xmax = 9.
t = 1.0
xi = np.linspace(xmin, xmax, num=nx)/t


# Load the C-library compiled with "make"
libc = cdll.LoadLibrary("./libstvenant_c.so")
# Explicit the C-function argument types
libc.riemann.argtypes = [np.ctypeslib.ndpointer(float, ndim=1,
                                                flags='aligned'),
                         np.ctypeslib.ndpointer(float, ndim=1,
                                                flags='aligned'),
                         c_double,
                         np.ctypeslib.ndpointer(float, ndim=1,
                                                flags='aligned, contiguous,'
                                                      'writeable')]


def riemann(wL, wR, xi):
    """A wrapper to the C-function libc.riemann"""

    wL = np.require(wL, float, ['ALIGNED'])
    wR = np.require(wR, float, ['ALIGNED'])
    w = np.zeros_like(wL)
    libc.riemann(wL, wR, xi, w)
    return w


def stvenant(plot_file=False):
    """main function that loops over x and plots the results"""
    start_time = default_timer()  # Start timer
    w = np.array([riemann(wL, wR, xi_j) for xi_j in xi])  # Loop over xi
    elapsed_time = default_timer() - start_time
    print("Time to run Riemann solver [s]: {}".format(elapsed_time))

    plt.plot(xi, w[:, 0], label="h")
    plt.plot(xi, w[:, 1]/w[:, 0], label="u")

    if plot_file:
        filename = "plotriem"
        data = pd.read_csv("plotriem", delim_whitespace=True,
                           header=None).values
        xi_f, h_f, u_f = data[:, 0], data[:, 1], data[:, 2]
        plt.plot(xi_f, h_f, 'g+', label="h_"+filename)
        plt.plot(xi_f, u_f, 'k+', label="u_"+filename)

    plt.legend()
    plt.show()


if __name__ == '__main__':
    stvenant(plot_file=True)
