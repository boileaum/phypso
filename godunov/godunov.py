#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
import godunov as godunov_python
from importlib import import_module
import matplotlib.pyplot as plt
import numpy as np
import timeit

PROBLEMS = set(('burgers', 'stvenant'))
KERNELS = set(('python', 'pythran', 'numpy', 'numba', 'fortran'))
CFL = 0.8
XMIN = -1.
XMAX = 2.


def burgers_exact(x, t):
    """return the exact solution at (x, t) of Burger's equation"""
    if t <= 1.:
        if x <= t: w = 1.
        elif x >= 1.: w = 0.
        elif (x > t) and (x <= 1.):
            w = (1. - x)/(1. - t)
    else:
        w = 1. if (x - 1.) <= 0.5*(t - 1.) else 0.
    return w


def L2_err(sol, ref, norm=1.0):
    """return normalized infinite-norm error between sol and ref"""
    N = len(sol)
    return np.linalg.norm(sol - ref, ord=2)/np.sqrt(N)/norm


class Problem():

    def __init__(self, **kwargs):

        self.problem = kwargs.get('problem', 'burgers')
        self.tmax = kwargs.get('tmax', 1.0)
        self.nmax = kwargs.get('nmax', 100)
        self.plot = kwargs.get('plot', False)
        self.profile = kwargs.get('profile', False)
        self.kernel = kwargs.get('kernel', 'python')

        print("problem:", self.problem)
        print("tmax =", self.tmax)
        print("nmax =", self.nmax)
        print("kernel:", self.kernel)

        if self.kernel == 'fortran':
            self.xm = np.zeros(self.nmax+2, order='F')
            self.wn = np.zeros(self.nmax+2, order='F')
        else:
            self.xm = np.zeros(self.nmax+2)
            self.wn = np.zeros(self.nmax+2)

        self.dx = float(XMAX - XMIN)/self.nmax
        self.cfl = CFL

        module_name = "godunov_" + self.kernel
        godunov_module = import_module(module_name)
        if self.problem == 'burgers':
            self.sol_exact = burgers_exact
            self.solver = godunov_module.Burgers()
        elif self.problem == 'stvenant':
            pass
        else:
            exit("Unknow problem:", self.problem)

        # Initialize with analytical solution
        self.xm = np.linspace(XMIN - 0.5*self.dx, XMAX + 0.5*self.dx,
                              num=self.nmax+2)
        self.wn = np.array([self.sol_exact(x, 0.) for x in self.xm])

    def compute_sol_num(self):
        self.wn = self.solver.timeloop(self.tmax, self.nmax,
                                       self.xm, self.wn, self.cfl)
#        self.exec_time = time

    def compute_sol_exact(self):
        self.wexact = np.vectorize(self.sol_exact)(self.xm, self.tmax)

    def compute_error(self):
        self.error = L2_err(self.wn, self.wexact)
        print("L2 error = {:f}".format(self.error))

    def solve(self):
        self.compute_sol_num()
        self.compute_sol_exact()
        self.compute_error()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solve hyperbolic problem")
    parser.add_argument('--problem', choices=PROBLEMS,
                        default='burgers', help="select Problem type")
    parser.add_argument('--tmax', type=float, metavar="final_time", default=1.,
                        help="simulation final time")
    parser.add_argument('--nmax', type=int, metavar="number_of_pts",
                        default=100, help="number of grid points")
    parser.add_argument('--profile', action='store_true',
                        help="activate profiling")
    parser.add_argument('--plot', action='store_true',
                        help="activate plotting")
    parser.add_argument('--kernel', choices=KERNELS,
                        default='python', help="select kernel type")
    args = parser.parse_args()

    p1 = Problem(**vars(args))
    p1.solve()

    if args.plot:
        plt.plot(p1.xm[1:-1], p1.wexact[1:-1], 'r+', label='exact')
        plt.plot(p1.xm[1:-1], p1.wn[1:-1], 'k-', label=p1.kernel)
        plt.legend()
        plt.show()
