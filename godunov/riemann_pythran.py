#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export riemann_burgers(float, float, float)
"""
Riemann solver
"""


def riemann_burgers(wL, wR, xi):
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
