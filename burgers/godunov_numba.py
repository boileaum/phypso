#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver for the Burgers equation
"""

from numba import jit
import numpy as np

xmin = -1.
xmax = 2.


def sol_exact(x, t):
    """return the exact solution at (x, t)"""
    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
    return w


@jit
def riemann(wl, wr, xi):
    """Return the solution of the Riemann problem at xi"""
    if wl > wr:
        sigma = 0.5*(wl + wr)
        if xi < sigma: w = wl
        if xi >= sigma: w = wr
    else:
        if xi <= wl: w = wl
        if xi >= wr: w = wr
        if (xi > wl) and (xi < wr): w = xi
    return w


def timeloop(tmax, nmax):
    """Iterate overt time to return the solution at t = tmax"""

    xm = np.zeros(nmax+2)
    wn = np.zeros(nmax+2)
    dx = float(xmax - xmin)/nmax
    cfl = 0.8

    # Initialization
    for i in range(nmax+2):
        xm[i] = xmin + (i - 0.5)*dx
        wn[i] = sol_exact(xm[i], 0.)

    wnp1 = wn.copy()

    t = 0.
    while t < tmax:
        vmax = max(abs(wn))
        dt = cfl*dx/vmax
        for i in range(1, nmax+1):
            # Right flux
            w = riemann(wn[i], wn[i+1], 0.)
            wnp1[i] -= 0.5*dt/dx*w**2
            # Left flux
            w = riemann(wn[i-1], wn[i], 0.)
            wnp1[i] += 0.5*dt/dx*w**2
        wn[:] = wnp1[:]  # Update
        t += dt

    return xm, wn
