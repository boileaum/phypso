#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export phys_flux(float[])
#pythran export numflux(float[], float[])
"""
Created on Tue Nov 14 00:17:12 2017
"""
import numpy
import riemann_stvenant

g = 9.81


def phys_flux(w):
    h = w[0]
    u = w[1]/h
    return numpy.array([h*u, h*u**2 + g*h**2/2.])


def numflux(wL, wR):
    """Return numerical flux using Riemann solver"""

    nmax = len(wL) - 1
    flux = numpy.zeros((nmax+1, 2))
    for i in range(nmax+1):
        w = riemann_stvenant.riemann(wL[i, :], wR[i, :], 0.)
        flux[i] = phys_flux(w)
    return flux
