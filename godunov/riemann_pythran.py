#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export riemann_burgers(float, float, float)
#pythran export numflux(float[], float[])
"""
Riemann solver
"""

import numpy as np


def riemann_burgers(wL, wR, xi):
    """Return the solution of the Riemann problem at xi"""
    if wL > wR:
        sigma = 0.5*(wL + wR)
        if xi < sigma: w = wL
        if xi >= sigma: w = wR
    else:
        if xi <= wL: w = wL
        if xi >= wR: w = wR
        if (xi > wL) and (xi < wR): w = xi
    return w


def numflux(wL, wR):
    """Return numerical flux using Riemann solver"""

    nmax = len(wL) - 1
    flux = np.zeros(nmax+1)
    for i in range(nmax+1):
        w = riemann_burgers(wL[i], wR[i], 0.)
        flux[i] = w*w/2.
    return flux
