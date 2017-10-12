#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
import godunov
import numpy as np
import timeit


def compute_sol(tmax, nmax):
    """Allocate arrays and solve Godunov problem"""
    xm = np.zeros(nmax+2)
    wn = np.zeros(nmax+2)
    wnp1 = np.zeros(nmax+2)
    xm, wn = godunov.timeloop(tmax, xm, wn, wnp1)
    return xm, wn


def L2_err(sim, ref, norm=1.0):
    """return normalized infinite-norm error between sim and ref"""
    N = len(sim)
    return np.linalg.norm(sim - ref, ord=2)/np.sqrt(N)/norm


def main(tmax, nmax, profile=False, plot=False):

    print("tmax =", tmax)
    print("nmax =", nmax)
    xm, wn = compute_sol(tmax, nmax)
    if profile:
        s = "xm, wn = compute_sol({}, {})".format(tmax, nmax)
        ntime = 10
        total_time = timeit.timeit(s, setup="import godunov", number=ntime,
                                   globals=globals())
        print("Mean time [s] over {} executions = {}".format(ntime,
              total_time/ntime))

    wex = godunov.sol_exact(tmax, xm)

    print("L2 error = {}".format(L2_err(wn, wex)))

#    with open('exact.dat', 'w') as f_exact:
#        with open('godu.dat', 'w') as f_godu:
#            for i in range(1, nmax+1):
#                f_exact.write("{} {}\n".format(xm[i], wex[i]))
#                f_godu.write("{} {}\n".format(xm[i], wn[i]))

    if plot:
        import matplotlib.pyplot as plt
        plt.plot(xm[1:-1], wex[1:-1], 'r+', label='exact')
        plt.plot(xm[1:-1], wn[1:-1], 'k-', label='godunov')
        plt.legend()


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
    args = parser.parse_args()

    main(args.tmax, args.nmax, args.profile, args.plot)
