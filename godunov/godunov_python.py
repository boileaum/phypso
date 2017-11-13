#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver
"""

import numpy as np


class Hyperbolic():
    """A generic class to define a hyperbolic solver"""

    def riemann(self):
        """Return the solution of the Riemann problem at xi"""
        pass

    def vmax(self, wn):
        pass

    def numflux(self):
        """Return numerical flux using Riemann solver"""
        pass

    def timeloop(self):
        """Iterate overt time to return the solution at t = tmax"""

        t = 0.
        while t < self.tmax:
            dt = self.cfl*self.dx/self.vmax()

            flux = self.numflux(self.wn[:-1], self.wn[1:])
            self.wn[1:-1] -= dt/self.dx*(flux[1:] - flux[:-1])

            t += dt

        return self.wn


class Burgers(Hyperbolic):

    cfl = 0.8
    xmin = -1.
    xmax = 2.

    def sol_exact(self, x, t):
        """return the exact solution at (x, t) of Burger's equation"""
        if t <= 1.:
            if x <= t: w = 1.
            elif x >= 1.: w = 0.
            elif (x > t) and (x <= 1.):
                w = (1. - x)/(1. - t)
        else:
            w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
        return w

    def init_arrays(self):
        self.xm = np.zeros(self.nmax+2)
        self.wn = np.zeros(self.nmax+2)

    def __init__(self, nmax, tmax):

        self.nmax = nmax
        self.tmax = tmax
        self.init_arrays()

        # Initialize with analytical solution
        self.dx = float(self.xmax - self.xmin)/self.nmax
        self.xm = np.linspace(self.xmin - 0.5*self.dx, self.xmax + 0.5*self.dx,
                              num=self.nmax+2)
        self.wn = np.array([self.sol_exact(x, 0.) for x in self.xm])

    def riemann(self, wL, wR, xi):
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

    def vmax(self):
        return max(abs(self.wn))

    def numflux(self, wL, wR):
        """Return numerical flux using Riemann solver"""

        nmax = len(wL) - 1
        flux = np.zeros(nmax+1)
        for i in range(nmax+1):
            w = self.riemann(wL[i], wR[i], 0.)
            flux[i] = w*w/2.
        return flux

    def compute_sol_exact(self):
        return np.vectorize(self.sol_exact)(self.xm, self.tmax)


