#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
import godunov as godunov_python
import numpy as np
import timeit


def compute_sol(tmax, nmax, kernel='python'):
    """Solve Godunov problem using the selected kernel"""
    if kernel == 'python':  # Use native (and naive!) python
        return godunov_python.timeloop(tmax, nmax)
    elif kernel == 'pythran':  # Use C++ with pythran
        import godunov_pythran
        return godunov_pythran.timeloop(tmax, nmax)
    elif kernel == 'numpy':
        import godunov_numpy
        return godunov_numpy.timeloop(tmax, nmax)
    elif kernel == 'fortran':
        import godunov_fortran
        # Allocate Fortran-ordered numpy arrays
        xm = np.zeros(nmax+2, order='F')
        wn = np.zeros(nmax+2, order='F')
        godunov_fortran.timeloop(tmax, xm, wn)
        return xm, wn


def L2_err(sol, ref, norm=1.0):
    """return normalized infinite-norm error between sol and ref"""
    N = len(sol)
    return np.linalg.norm(sol - ref, ord=2)/np.sqrt(N)/norm


def main(tmax, nmax, profile=False, plot=False, kernel='python'):

    print("tmax =", tmax)
    print("nmax =", nmax)
    print("kernel:", kernel)

    xm, wn = compute_sol(tmax, nmax, kernel)  # Always compute once
    if profile:
        s = "xm, wn = compute_sol({}, {}, kernel='{}')".format(tmax, nmax,
                                                               kernel)
        ntime = 10
        total_time = timeit.timeit(s, number=ntime, globals=globals())

        print("Mean time [s] over {} executions = {}".format(ntime,
              total_time/ntime))

    wex = np.vectorize(godunov_python.sol_exact)(xm, tmax)

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
                        default=100, help="number of grid points")
    parser.add_argument('--profile', action='store_true',
                        help="activate profiling")
    parser.add_argument('--plot', action='store_true',
                        help="activate plotting")
    parser.add_argument('--kernel', choices=['python', 'pythran', 'numpy',
                                             'fortran'],
                        default='python', help="select kernel type")
    args = parser.parse_args()

    main(args.tmax, args.nmax, args.profile, args.plot, args.kernel)
