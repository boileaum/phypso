#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export sol_exact(float, float[])
#pythran export riemann(float, float, float)
#pythran export timeloop(float, float[], float[], float[])
#pythran export A1, A2
"""
Created on Wed Oct  4 17:10:53 2017

@author: boileau
"""

import numpy as np

A1 = -1.
A2 = 2.

def sol_exact_x(x, t):

    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
    return w


def sol_exact(tmax, xm):

    wex = np.zeros_like(xm)
    for i in range(1, len(wex)-1):
        wex[i] = sol_exact_x(xm[i], tmax)
    return wex


def riemann(wl, wr, xi):

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

    nmax = len(xm) - 2
    dx = float(A2 - A1)/nmax
    cfl = 0.8

    # Initialization
    for i in range(nmax+1):
        xm[i] = A1 + (i + 0.5)*dx
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
        print("t =", t)

    return xm, wn
