#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver, numpy version
"""

from godunov_python import Burgers as Burgers_python


class Burgers(Burgers_python):
    """A set of methods to solve Burgers' problem using a numpy kernel"""

    def riemann(self, wL, wR, xi):
        """Exact Riemann solver for Burgers' equation"""
        sigma = (wL + wR)/2.
        w = (wL < wR)*((xi < wL)*wL +
                       (xi > wR)*wR + ((xi >= wL) & (xi <= wR))*xi) + \
            (wL >= wR)*((xi < sigma)*wL + (xi >= sigma)*wR)
        return w

    def numflux(self, wL, wR):
        """Return numerical flux using Riemann solver"""
        w = self.riemann(wL, wR, 0.)
        flux = w*w/2.
        return flux
