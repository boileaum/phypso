#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver for the Burgers equation using numpy
(This numpy function is in a separate module because pythran does not fully
support numpy slicing)
"""

from godunov import xmin, xmax, sol_exact
import numpy as np


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

        wl = wn[:-1]
        wr = wn[1:]

        sigma = 0.5*(wl + wr)
        w_shock = wl if np.where(sigma > 0.) else wr

        condlist = [wl >= 0., wr <= 0., (wl < 0) * (wr > 0)]
        choicelist = [wl, wr, 0.]
        w_expan = np.select(condlist, choicelist)

        w = w_shock if np.where(wl > wr) else w_expan

        wn[1:-1] = wn[1:-1] - 0.5*dt/dx*w[1:]**2 + 0.5*dt/dx*w[:-1]**2
        t += dt

    return xm, wn
