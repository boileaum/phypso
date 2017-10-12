#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export sol_exact_x(float, float)
#pythran export sol_exact(float, float[])
#pythran export riemann(float, float, float)
#pythran export timeloop(float, float[], float[], float[])
#pythran export xmin, xmax
"""
Godunov solver for the Burgers equation
"""

import numpy as np

xmin = -1.
xmax = 2.


def sol_exact_x(x, t):
    """return exact solution for at (x, t)"""
    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
    return w


def sol_exact(tmax, xm):
    """loop on each cell to return the exact solution"""
    wex = np.zeros_like(xm)
    for i in range(1, len(wex)-1):
        wex[i] = sol_exact_x(xm[i], tmax)
    return wex


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


def timeloop(tmax, xm, wn, wnp1):
    """Iterate overt time to return the solution at t = tmax"""
    nmax = len(xm) - 2
    dx = float(xmax - xmin)/nmax
    cfl = 0.8

    # Initialization
    for i in range(nmax+1):
        xm[i] = xmin + (i + 0.5)*dx
        wn[i] = sol_exact_x(xm[i], 0.)
        wnp1[i] = wn[i]

    t = 0.
    while t < tmax:
        vmax = max(abs(wn))
        dt = cfl*dx/(vmax + 1e-16)
        for i in range(1, nmax+1):
            # Right flux
            w = riemann(wn[i], wn[i+1], 0.)
            wnp1[i] = wnp1[i] - 0.5*dt/dx*w**2
            # Left flux
            w = riemann(wn[i-1], wn[i], 0.)
            wnp1[i] = wnp1[i] + 0.5*dt/dx*w**2
        wn[:] = wnp1[:]  # Update
        t += dt
        #print("t =", t)

    return xm, wn
