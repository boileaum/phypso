#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export sol_exact(float, float)
#pythran export riemann(float, float, float)
#pythran export timeloop(float, int)
#pythran export xmin, xmax
"""
Godunov solver for the Burgers equation
"""

import numpy as np
#from TP_EDP3.phypso.stvenant.riemann import riemann as riemann_stvenant
from riemann_stvenant import riemann as riemann_stvenant
import matplotlib.pyplot as plt
# xmin = -1.
# xmax = 2.

g = 9.81

m = 2
nx = 1001

xmin = -10.
xmax = 10.
t = 1.0
#xi = np.linspace(xmin, xmax, num=nx)/t




def sol_exact(x, t):
    """return the exact solution at (x, t)"""
    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
    return w


def sol_exact_stvenant(x,t):
    t += 1e-10
    xi=x/t
    hL, uL = 2., 0.
    hR, uR = 1., 0.

    # wL = np.array([hL, hL * uL])
    # wR = np.array([hR, hR * uR])
    wL = np.array([1., -0.5]) #pour le deuxième cas
    wR = np.array([1., 0.5])

    return riemann_stvenant(wL,wR,xi)


# noinspection PyUnboundLocalVariable
def riemann(wl, wr, xi):
    """Return the solution of the Riemann problem at xi"""
    if wl > wr:
        sigma = 0.5*(wl + wr)
        if xi < sigma: w = wl
        if xi >= sigma: w = wr
    else:
        if xi <= wl: w = wl
        if xi >= wr: w = wr
        if (xi > wl) and (xi < wr): w = xi
    return w

def vmax(w):
    return abs(w)

def vmax_stvenant(w):
    h = w[0]
    u = w[1]/h
    return np.sqrt(g*h) +np.abs(u)

def flux_phy(w):
    return w**2/2


def flux_phy_stvenant(w):
    h = w[0]
    u =w[1]/h
    return np.array([h*u,h*u*u+g*h*h/2])


def timeloop(tmax, nmax,riemann=riemann,sol_exact=sol_exact,vmax=vmax
             ,flux_phy=flux_phy,plot=False,wall=False):
    """Iterate overt time to return the solution at t = tmax"""

    xm = np.zeros(nmax+2)
    wn = np.zeros((nmax+2,m))
    dx = float(xmax - xmin)/nmax
    cfl = 0.8

    # Initialization
    for i in range(nmax+2):
        xm[i] = xmin + (i - 0.5)*dx
        wn[i,:] = sol_exact(xm[i], 0.)

    wnp1 = wn.copy()

    t = 0.
    while t < tmax:
        v = np.max([vmax(wn[i]) for i in range(nmax+2)])
        #vmax = np.apply_along_axis()
        dt = cfl*dx/v
        #print("t=",t,"dt=",dt,"v=",v)
        mass = np.sum(wn[1:-1,0]*dx)
        print("mass= ",mass)#pour vérifier la conservation de la mass
        if wall:
            wn[0,0] = wn[1,0]
            wn[0,1] = -wn[1,1]
            wn[nmax+1,0] = wn[nmax,0]
            wn[nmax+1,1] = -wn[nmax,1]
        for i in range(1, nmax+1):
            # Right flux
            w = riemann(wn[i], wn[i+1], 0.)
            flux = flux_phy(w)
            wnp1[i] -= dt/dx*flux
            # Left flux
            w = riemann(wn[i-1], wn[i], 0.)
            flux = flux_phy(w)
            wnp1[i] += dt/dx*flux

        wn[:] = wnp1[:]  # Update
        t += dt
        if plot:
            plt.clf()
            plt.plot(xm[1:-1], wn[1:-1,0])
            plt.ylim(0,2.5)
            plt.draw()
            plt.pause(0.001)

    wnexact=np.zeros((nmax+2,m))
    """pour obtenir la valeur exacte à la fin"""
    for i in range(nmax+2):
        wnexact[i,:]=sol_exact(xm[i],t)
    return xm, wn ,wnexact

tmax = 10./4.4
nmax = 200
xm, wn , wnexact = timeloop(tmax,nmax
                  ,riemann=riemann_stvenant
                  ,sol_exact=sol_exact_stvenant
                  ,vmax=vmax_stvenant
                    ,flux_phy=flux_phy_stvenant
                    ,plot=False
                            ,wall=True)
plt.clf()
"""affichage de la valeur approché et de la valeur exacte"""
# plt.plot(xm[1:-1],wn[1:-1,0],label="approx")
# plt.plot(xm[1:-1],wnexact[1:-1,0],label="exacte")

"""affichage de la hauteur et de la vitesse"""
plt.plot(xm[1:-1],wnexact[1:-1,0],label="hauteur de l'eau")
plt.plot(xm[1:-1],wnexact[1:-1,1]/wnexact[1:-1,0],label="vitesse")
plt.legend()
plt.show()
