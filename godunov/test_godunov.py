#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that numerical error of the Godunov solver with various kernel versions
"""

from godunov import Problem, PROBLEMS
from pytest import mark, approx
import itertools


nmax_kernel_combinations = itertools.product([100, 1000], PROBLEMS['burgers'])


@mark.parametrize('nmax, kernel', nmax_kernel_combinations)
def test_burgers(nmax, kernel):
    """Test a combination of nmax values and kernel versions and compare to
    expected numerical error"""
    tmax = 1.0
    err_ref = {100: 0.04496958454369648,
               1000: 0.018877105158683641}
    problem = Problem(tmax=tmax, nmax=nmax, kernel=kernel, problem="burgers")
    problem.solve()
    assert problem.error == approx(err_ref[nmax])


nmax_kernel_combinations = itertools.product([100, 500], PROBLEMS['stvenant'])


@mark.parametrize('nmax, kernel', nmax_kernel_combinations)
def test_stvenant(nmax, kernel):
    """Test a combination of nmax values and kernel versions and compare to
    expected numerical error"""
    tmax = 1.0
    err_ref = {100: 0.124886436882,
               500: 0.0562642803228}
    problem = Problem(tmax=tmax, nmax=nmax, kernel=kernel, problem="stvenant")
    problem.solve()
    assert problem.error == approx(err_ref[nmax])
