#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve St-Venant's problem in 1D calling a Riemann solver in C
(See https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.python-as-glue.html#index-3)
"""

from ctypes import cdll, c_double
import matplotlib.pyplot as plt
import numpy as np


hL = 2.
uL = 0.
hR = 1.
uR = 0.

wL = np.array([hL, hL*uL])
wR = np.array([hR, hR*uR])

nx = 1000

xmin = -9.
xmax = 9.
dx = (xmax - xmin)/nx

w = np.zeros((nx, 2))
wtype = w[0].dtype
x = np.linspace(xmin, xmax, num=nx)

t = 1.

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


def stvenant():

    w = np.array([riemann(wL, wR, xi) for xi in x])  # Loop over xi

    plt.plot(x, w[:, 0], label="h")
    plt.plot(x, w[:, 1]/w[:, 0], label="u")
    plt.legend()
    plt.show()


if __name__ == '__main__':
    stvenant()
