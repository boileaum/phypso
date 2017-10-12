#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
import godunov
import numpy as np
import timeit


def compute_sol_python(tmax, nmax):
    """Allocate arrays and solve Godunov problem using Python (or C++ through
    Pythran)"""
    xm = np.zeros(nmax+2)
    wn = np.zeros(nmax+2)
    wnp1 = np.zeros(nmax+2)
    xm, wn = godunov.timeloop(tmax, xm, wn, wnp1)
    return xm, wn


def compute_sol_fortran(tmax, nmax):
    """Allocate Fortran-ordered arrays and solve Godunov problem using
    Fortran"""
    import libgodunov
    xm = np.zeros(nmax+2, order='F')
    wn = np.zeros(nmax+2, order='F')
    libgodunov.timeloop(tmax, xm, wn)
    return xm, wn


def L2_err(sol, ref, norm=1.0):
    """return normalized infinite-norm error between sol and ref"""
    N = len(sol)
    return np.linalg.norm(sol - ref, ord=2)/np.sqrt(N)/norm


def main(tmax, nmax, profile=False, plot=False, fortran=False):

    print("tmax =", tmax)
    print("nmax =", nmax)
    global compute_sol
    compute_sol = compute_sol_fortran if fortran else compute_sol_python
    xm, wn = compute_sol(tmax, nmax)  # Always compute once
    if profile:
        s = "xm, wn = compute_sol({}, {})".format(tmax, nmax)
        ntime = 10
        total_time = timeit.timeit(s, setup="import godunov", number=ntime,
                                   globals=globals())
        print("Mean time [s] over {} executions = {}".format(ntime,
              total_time/ntime))

    wex = godunov.sol_exact(tmax, xm)

    print("L2 error = {}".format(L2_err(wn, wex)))

    if plot:
        import matplotlib.pyplot as plt
        plt.plot(xm[1:-1], wex[1:-1], 'r+', label='exact')
        plt.plot(xm[1:-1], wn[1:-1], 'k-', label='godunov')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solve Burgers problem")
    parser.add_argument('--tmax', type=float, metavar="final_time", default=1.,
                        help="simulation final time")
    parser.add_argument('--nmax', type=int, metavar="number_of_pts",
                        default=100, help="number of points")
    parser.add_argument('--profile', action='store_true',
                        help="activate profiling")
    parser.add_argument('--plot', action='store_true',
                        help="activate plotting")
    parser.add_argument('--fortran', action='store_true',
                        help="use fortran binding")
    args = parser.parse_args()

    main(args.tmax, args.nmax, args.profile, args.plot, args.fortran)
