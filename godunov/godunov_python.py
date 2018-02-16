#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver
"""

import numpy as np
from burgers.riemann import riemann as riemann_burgers
from stvenant.riemann import riemann as riemann_stvenant
from stvenant.riemann import g


class Hyperbolic():
    """A generic class to define a hyperbolic solver"""

    def __init__(self, problem):

        self.problem = problem
        self.nmax = self.problem.nmax
        self.tmax = self.problem.tmax
        self.plot = self.problem.plot

        # Initialize with analytical solution
        self.dx = float(self.xmax - self.xmin)/self.nmax
        self.xm = np.linspace(self.xmin - 0.5*self.dx, self.xmax + 0.5*self.dx,
                              num=self.nmax+2)
        self.init_sol()
        self.init_mass = self.problem.get_mass(self.wn, self.dx)

    def sol_exact(self, x, t):
        pass

    def init_sol(self):
        self.wn = np.array([self.sol_exact(x, 0.) for x in self.xm])

    def riemann(self):
        """Return the solution of the Riemann problem at xi"""
        pass

    def vmax(self, wn):
        pass

    def phys_flux(self, w):
        pass

    def numflux(self, wL, wR):
        """Return numerical flux using Riemann solver"""

        for i in range(self.nmax+1):
            w = self.riemann(wL[i], wR[i], 0.)
            self.flux[i] = self.phys_flux(w)
        return self.flux

    def compute_sol_exact(self, t):
        return np.array([self.sol_exact(xm, t) for xm in self.xm])

    def timeloop(self):
        """Iterate overt time to return the solution at t = tmax"""

        if self.plot:
            import matplotlib.pyplot as plt
            plt.ion()
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_title(self.problem.message)
            self.problem.plot_wn(self.ax, self.wn)

        t = 0.
        while t < self.tmax:
            self.problem.BC(self.wn)
            dt = self.cfl*self.dx/self.vmax()

            self.flux = self.numflux(self.wn[:-1], self.wn[1:])
            self.wn[1:-1] -= dt/self.dx*(self.flux[1:] - self.flux[:-1])

            t += dt
            if self.plot:
                self.problem.plot_update(self.ax, self.wn)
                plt.show()
                plt.pause(0.001)
            self.mass = self.problem.get_mass(self.wn, self.dx)
            # print("Total mass =", self.mass)

        if self.plot:
            plt.ioff()
        return self.wn


class Burgers(Hyperbolic):

    cfl = 0.8
    xmin = -1.
    xmax = 2.

    def __init__(self, problem):

        self.riemann = riemann_burgers
        super().__init__(problem)
        self.flux = np.zeros(self.nmax+1)

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

    def vmax(self):
        return max(abs(self.wn))

    def phys_flux(self, w):
        return w**2/2.


class StVenant(Hyperbolic):

    cfl = 0.8
    xmin = -10.
    xmax = 10.
    hL, uL = 2., 0.
    hR, uR = 1., 0.
    wL = np.array([hL, hL*uL])
    wR = np.array([hR, hR*uR])

    def __init__(self, problem):

        self.riemann = riemann_stvenant
        super().__init__(problem)
        self.flux = np.zeros((self.nmax+1, 2))

    def sol_exact(self, x, t):
        """return the exact solution at (x, t) of Burger's equation"""
        t += 1.e-10
        xi = x/t
        w = self.riemann(self.wL, self.wR, xi)
        return w

    def vel(self, w):
        h = w[0]
        u = w[1]/h
        return np.sqrt(g*h) + np.abs(u)

    def vmax(self):
        return np.max([self.vel(self.wn[i, :]) for i in range(self.nmax)])

    def phys_flux(self, w):
        h = w[0]
        u = w[1]/h
        return np.array([h*u, h*u**2 + g*h**2/2.])
