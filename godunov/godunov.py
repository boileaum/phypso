#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
from importlib import import_module
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer

PROBLEMS = set(('burgers', 'stvenant'))
KERNELS = set(('python', 'numpy', 'pythran'))


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

        module_name = "godunov_" + self.kernel
        godunov_module = import_module(module_name)
        if self.problem == 'burgers':
            self.solver = godunov_module.Burgers(self.nmax, self.tmax)
        elif self.problem == 'stvenant':
            pass
        else:
            exit("Unknow problem:", self.problem)

        self.xm = self.solver.xm

    def __repr__(self):

        return "problem: {}\n".format(self.problem) + \
               "tmax = {}\n".format(self.problem) + \
               "nmax = {}\n".format(self.problem) + \
               "kernel: {}".format(self.kernel)

    def compute_sol_num(self):

        start_time = default_timer()  # Start timer
        self.wn = self.solver.timeloop()
        self.exec_time = default_timer() - start_time
        print("Time to solve {} with {} [s]: {:f}".format(self.problem,
                                                          self.kernel,
                                                          self.exec_time))

    def compute_sol_exact(self):
        self.wexact = self.solver.compute_sol_exact()

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

    p = Problem(**vars(args))
    print(p)
    p.solve()

    if args.plot:
        plt.plot(p.xm[1:-1], p.wexact[1:-1], 'r+', label='exact')
        plt.plot(p.xm[1:-1], p.wn[1:-1], 'k-', label=p.kernel)
        plt.legend()
        plt.show()
