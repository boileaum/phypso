#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 17:10:53 2017

@author: boileau
"""

import numpy as np

NMAX = 100
TMAX = 1.5
A1 = -1.
A2 = 2.
xm = np.zeros(NMAX+2)
wex = np.zeros(NMAX+2)
wn = np.zeros(NMAX+2)
wnp1 = np.zeros(NMAX+2)


def sol_exact_x(x, t):

    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.

    return w


def sol_exact():
    for i in range(1, NMAX+1):
        wex[i] = sol_exact_x(xm[i], TMAX)
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


def timeloop():

    dx = (A2 - A1)/NMAX
    cfl = 0.8

    # Initialization
    for i in range(NMAX+1):
        xm[i] = A1 + (i + 0.5)*dx
        wn[i] = sol_exact_x(xm[i], 0.)
        wnp1[i] = wn[i]

    t = 0.
    while t < TMAX:
        vmax = max(abs(wn))
        dt = cfl*dx/vmax
        for i in range(1, NMAX+1):
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
