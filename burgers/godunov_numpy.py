#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver for the Burgers equation using numpy
(This numpy function is in a separate module because pythran does not fully
support numpy slicing)
"""

from godunov import xmin, xmax, sol_exact
import numpy as np


def riemann_burgers(wL, wR, xi):
    """Exact Riemann solver for Burgers' equation"""
    sigma = (wL + wR)/2.
    w = (wL < wR)*((xi < wL)*wL +
                   (xi > wR)*wR + ((xi >= wL) & (xi <= wR))*xi) + \
        (wL >= wR)*((xi < sigma)*wL + (xi >= sigma)*wR)
    return w


def numflux(wL, wR):
    """Return numerical flux using Riemann solver"""
    w = riemann_burgers(wL, wR, 0.)
    flux = w*w/2.
    return flux


def timeloop(tmax, nmax):
    """Iterate overt time to return the solution at t = tmax"""
    dx = float(xmax - xmin)/nmax
    cfl = 0.8

    xm = np.linspace(xmin - 0.5*dx, xmax + 0.5*dx, num=nmax+2)
    wn = np.fromiter((sol_exact(x, 0.) for x in xm), xm.dtype)

    t = 0.
    while t < tmax:
        vmax = max(abs(wn))
        dt = cfl*dx/vmax

        flux = numflux(wn[:-1], wn[1:])
        wn[1:-1] -= dt/dx*(flux[1:] - flux[:-1])

        t += dt

    return xm, wn
