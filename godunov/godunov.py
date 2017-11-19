#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve Burgers' equation using the 1rst-order Godunov method
"""

import argparse
from importlib import import_module
import matplotlib.pyplot as plt
import numpy as np
import pandas
from timeit import default_timer

all_KERNELS = 'python', 'numpy', 'pythran'
PROBLEMS = 'burgers', 'stvenant'


def L2_err(sol, ref, norm=1.0):
    """return normalized infinity norm error between sol and ref"""
    N = len(sol)
    return np.linalg.norm(sol - ref, ord=2)/np.sqrt(N)/norm


class Burgers():

    problem = "burgers"
    kernels = all_KERNELS

    def set_solver(self):
        self.solver = self.godunov_module.Burgers(self)

    def __init__(self, **kwargs):

        self.tmax = kwargs.get('tmax', 1.0)
        self.nmax = kwargs.get('nmax', 100)
        self.plot = kwargs.get('plot', False)
        self.profile = kwargs.get('profile', False)
        self.kernel = kwargs.get('kernel', 'python')

        # Import module according to kernel type
        module_name = "godunov_" + self.kernel
        self.godunov_module = import_module(module_name)
        self.set_solver()

    def __repr__(self):

        return "problem: {}\n".format(self.problem) + \
               "tmax = {}\n".format(self.tmax) + \
               "nmax = {}\n".format(self.nmax) + \
               "kernel: {}".format(self.kernel)

    def compute_sol_num(self):

        start_time = default_timer()  # Start timer
        self.wn = self.solver.timeloop()
        self.exec_time = default_timer() - start_time
        print("Time to solve {} with {} [s]: {:f}".format(self.problem,
                                                          self.kernel,
                                                          self.exec_time))

    def compute_sol_exact(self):
        self.wexact = self.solver.compute_sol_exact(self.tmax)

    def compute_error(self):
        self.error = L2_err(self.wn, self.wexact)
        print("L2 error = {:f}".format(self.error))

    def solve(self):
        self.compute_sol_num()
        self.compute_sol_exact()
        self.compute_error()

    def plot_exact(self, ax):
        x = self.xm[1:-1]
        ax.plot(x, self.wexact[1:-1], 'r+', label='exact')

    def plot_num(self, ax):
        x = self.xm[1:-1]
        ax.plot(x, self.wn[1:-1], label='{}, {}'.format(self.nmax,
                self.kernel))


class StVenant(Burgers):

    problem = "stvenant"
    kernels = 'python', 'pythran'

    def set_solver(self):
        self.solver = self.godunov_module.StVenant(self)

    def plot_exact(self, ax):
        x = self.xm[1:-1]
        h = self.wn[1:-1, 0]
        u = self.wn[1:-1, 1]/h
        ax.plot(x, h, '+', label='h_exact')
        ax.plot(x, u, '+', label='u_exact')

    def plot_num(self, ax):
        x = self.xm[1:-1]
        h = self.wn[1:-1, 0]
        u = self.wn[1:-1, 1]/h
        ax.plot(x, h, label='h ({}, {})'.format(self.nmax, self.kernel))
        ax.plot(x, u, label='u ({}, {})'.format(self.nmax, self.kernel))


def plot(*problems):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot exact solution only for first problem
    p0 = problems[0]
    p0.plot_exact(ax)

    for p in problems:
        p.plot_num(ax)

    ax.set_xlabel(r'$x$')
    fig.legend()
    plt.show()


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
    parser.add_argument('--kernel', choices=all_KERNELS + ("bench",),
                        default='python', help="select kernel type")
    args = parser.parse_args()

    if args.problem == 'burgers':
        Problem = Burgers
    elif args.problem == 'stvenant':
        Problem = StVenant

    if args.kernel == "bench":
        del(args.kernel)
        del(args.nmax)
        sizes = 100, 500, 1000
        scores = pandas.DataFrame(data=0, columns=Problem.kernels,
                                  index=sizes)
        problems = []
        for nmax in sizes:
            for kernel in Problem.kernels:
                p = Problem(kernel=kernel, nmax=nmax, **vars(args))
                problems.append(p)
                p.solve()
                scores.loc[nmax, kernel] = p.exec_time

        if args.plot:
            plot(*problems)
        print("Execution times [s] for {} problem.".format(args.problem))
        print(scores)

        normalized_scores = scores.copy()
        for column in normalized_scores.columns:
            normalized_scores[column] /= scores['python']
        print("\nNormalized execution times [-]")
        print(normalized_scores)

    else:
        if args.kernel not in Problem.kernels:
            msg = "{} kernel not implemented for {} problem.".format(
                    args.kernel, args.problem)
            exit(msg)

        p = Problem(**vars(args))
        print(p)
        p.solve()

        if args.plot:
            plot(p)
