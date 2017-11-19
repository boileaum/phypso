#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export g
#pythran export Z(float, float)
#pythran export riemann(float[], float[], float)

"""
Riemann solver
"""

from math import sqrt
import numpy as np

g = 9.81


def Z(h1, h2):
    if h1 < h2:
        return 2.*sqrt(g)/(sqrt(h1) + sqrt(h2))
    else:
        return sqrt(g*(h1 + h2)/h1/h2/2.)


def riemann(wL, wR, xi):

    def fz(hL, hR, uL, uR, hs):

        if hs < hR:
            a = 0.6264183906e1/(sqrt(hs) + sqrt(hR))
        else:
            a = sqrt((4.905*hs + 4.905*hR)/hs/hR)

        if hs < hL:
            b = 0.6264183906e1/(sqrt(hs) + sqrt(hL))
        else:
            b = sqrt((4.905*hs + 4.905*hL)/hs/hL)

        return (hs - hR)*a + (hs - hL)*b - uL + uR

    def dfz(hL, hR, uL, uR, hs):
        if hs < hR:
            a = 0.6264183906e1/(sqrt(hs) + sqrt(hR))
        else:
            a = sqrt((4.905*hs + 4.905*hR)/hs/hR)

        if hs < hR:
            b = -0.3132091953e1*pow(sqrt(hs) + sqrt(hR), -0.2e1) \
                * pow(hs, -0.1e1 / 0.2e1)
        else:
            b = pow((4.905*hs + 4.905*hR)/hs/hR,
                    -0.1e1/0.2e1)*(4.905/hs/hR - (4.905*hs + 4.905*hR)
                                   * pow(hs, -0.2e1)/hR)/0.2e1

        if hs < hL:
            c = 0.6264183906e1/(sqrt(hs) + sqrt(hL))
        else:
            c = sqrt((4.905*hs + 4.905*hL) / hs / hL)

        if hs < hL:
            d = -0.3132091953e1*pow(sqrt(hs) + sqrt(hL), -0.2e1) \
                * pow(hs, -0.1e1/0.2e1)
        else:
            d = pow((4.905*hs + 4.905*hL)/hs/hL,
                    -0.1e1/0.2e1)*(4.905/hs/hL - (4.905*hs + 4.905*hL)
                                   * pow(hs, -0.2e1)/hL)/0.2e1

        return a + (hs - hR)*b + c + (hs - hL)*d

    hL = wL[0]
    hR = wR[0]
    uL = wL[1]/hL
    uR = wR[1]/hR

    hs = 1e-8
    dh = 1e8

#    crit = (uR - uL < 2*sqrt(g)*(sqrt(hL) - sqrt(hR)))
#    assert(crit)  # Void appears

    while(abs(dh) > 1e-6):

        f = fz(hL, hR, uL, uR, hs)
        df = dfz(hL, hR, uL, uR, hs)
        dh = -f / df
        hs += dh

    us = uR + (hs - hR)*Z(hs, hR)

    # left wave
    if hs < hL:  # 1-expansion
        v1m = uL - sqrt(g*hL)
        v1p = us - sqrt(g*hs)
    else:  # 1-chock
        hc = (hs + hL)/2
        alpha = sqrt(hs)/(sqrt(hs) + sqrt(hL))
        uc = alpha*us + (1 - alpha)*uL
        v1m = uc - sqrt(g*hc)
        v1p = v1m

    # right wave
    if hs < hR:  # 2-expansion
        v2m = us + sqrt(g*hs)
        v2p = uR + sqrt(g*hR)
    else:  # 2-shock
        hc = (hs + hR)/2
        alpha = sqrt(hs)/(sqrt(hs) + sqrt(hR))
        uc = alpha * us + (1 - alpha)*uR
        v2m = uc + sqrt(g*hc)
        v2p = v2m

    if xi < v1m:
        u = uL
        h = hL
    elif xi < v1p:
        u = (2*xi + uL + 2*sqrt(g*hL))/3
        h = (u - xi)*(u - xi)/g
    elif xi < v2m:
        u = us
        h = hs
    elif xi < v2p:
        u = (2*xi + uR - 2*sqrt(g*hR))/3
        h = (u - xi)*(u - xi)/g
    else:
        u = uR
        h = hR

    return np.array([h, h*u])
