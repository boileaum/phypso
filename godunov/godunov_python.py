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

    def timeloop(self, tmax, nmax, xm, wn, cfl):
        """Iterate overt time to return the solution at t = tmax"""

        dx = xm[1] - xm[0]
        t = 0.
        while t < tmax:
            dt = cfl*dx/self.vmax(wn)

            flux = self.numflux(wn[:-1], wn[1:])
            wn[1:-1] -= dt/dx*(flux[1:] - flux[:-1])

            t += dt

        return wn


class Burgers(Hyperbolic):

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

    def vmax(self, wn):
        return max(abs(wn))

    def numflux(self, wL, wR):
        """Return numerical flux using Riemann solver"""

        nmax = len(wL) - 1
        flux = np.zeros(nmax+1)
        for i in range(nmax+1):
            w = self.riemann(wL[i], wR[i], 0.)
            flux[i] = w*w/2.
        return flux
