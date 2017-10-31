#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve St-Venant's problem in 1D by calling a Riemann solver written in C.
See:
https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.python-as-glue.html#index-3
"""

import argparse
from ctypes import cdll, c_double
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from timeit import default_timer
import functools
import subprocess
import sys


KERNELS = ['python', 'C-lib', 'cython', 'pythran']


hL, uL = 2., 0.
hR, uR = 1., 0.

wL = np.array([hL, hL*uL])
wR = np.array([hR, hR*uR])

nx = 1001

xmin = -9.
xmax = 9.
t = 1.0
xi = np.linspace(xmin, xmax, num=nx)/t


def timer(func):
    """Decorator that display the function execution time"""
    @functools.wraps(func)
    def func_call(*args, **kwargs):
        start_time = default_timer()  # Start timer
        result = func(*args, **kwargs)
        elapsed_time = default_timer() - start_time
        print("Time to run {} with {} [s]: {:f}".format(func.__name__,
                                                        args[0],
                                                        elapsed_time))
        return result
    return func_call


# Load the C-library compiled with "make"

lib_filename = os.path.join(os.path.dirname(__file__), 'libstvenant_c.so')
libc = cdll.LoadLibrary(lib_filename)
# Explicit the C-function argument types
libc.riemann.argtypes = [np.ctypeslib.ndpointer(float, ndim=1,  # wL
                                                flags='aligned'),
                         np.ctypeslib.ndpointer(float, ndim=1,  # wR
                                                flags='aligned'),
                         c_double,  # xi_j
                         np.ctypeslib.ndpointer(float, ndim=1,  # w
                                                flags='aligned, contiguous,'
                                                      'writeable')]


def riemann_C(wL, wR, xi_j):
    """A wrapper to the C-function libc.riemann"""

    wL = np.require(wL, float, ['ALIGNED'])
    wR = np.require(wR, float, ['ALIGNED'])
    wi = np.zeros_like(wL)
    libc.riemann(wL, wR, xi_j, wi)
    return wi


@timer
def riemann_loop(kernel, wL, wR, xi):
    """Loop over xi to return w as a numpy array of size nx using the
    Riemann solver provided by the riemann_func function"""
    if kernel == "python":
        from riemann import riemann
    elif kernel == "cython":
        from cstvenant import riemann_cython as riemann
    elif kernel == "pythran":
        from riemann_pythran import riemann
    elif kernel == "C-lib":
        riemann = riemann_C
    else:
        sys.exit("Unknown kernel")
    return np.array([riemann(wL, wR, xi_j) for xi_j in xi])


def load_file(filename="plotriem"):
    """Reads file and return a (x, h, u) tuple"""
    file_path = os.path.join(os.path.dirname(__file__), filename)
    data = pd.read_csv(file_path, delim_whitespace=True, header=None).values
    return data[:, 0], data[:, 1], data[:, 2]


@timer
def C_program(name):
    """Run the full C program"""
    cwd = os.getcwd()
    program_dir = os.path.dirname(__file__)
    os.chdir(program_dir)

    cmd = [name] + ['{}'.format(arg) for arg in
                    (hL, uL, hR, uR, nx - 1, xmin, xmax, t)]

    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print("\033[31mError in C program {}".format(cmd[0]))
        print("Command line:\n"
              "  {}\n"
              "returns:".format(" ".join(cmd)))
        print(e.stderr.decode(), '\033[0m')

    os.chdir(cwd)


def stvenant(plot=False):
    """main function that loops over x and plots the results"""
    C_program("./stvenant.exe")

    if sys.platform == "darwin":
        print("Skip pythran kernel because it does not work on Mac currently")
        KERNELS.remove("pythran")

    for kernel in KERNELS:
        w = riemann_loop(kernel, wL, wR, xi)
        if plot:
            plt.plot(xi, w[:, 0], label="h_"+kernel)
            plt.plot(xi, w[:, 1]/w[:, 0], label="u_"+kernel)

    if plot:
        filename = "plotriem"
        xi_f, h_f, u_f = load_file(filename)
        plt.plot(xi_f, h_f, 'g+', label="h_"+filename)
        plt.plot(xi_f, u_f, 'k+', label="u_"+filename)

        plt.xlabel(r'$\xi$')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solve 1D St-Venant problem")
    parser.add_argument('--noplot', action='store_true',
                        help="Do not plot solution")
    args = parser.parse_args()
    stvenant(not args.noplot)
