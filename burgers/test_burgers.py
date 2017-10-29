#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that numerical error of the Godunov solver with various kernel versions
"""

from burgers import burgers, KERNELS
from pytest import mark, approx, param
import itertools
import sys


cartesian_product = itertools.product([100, 1000], KERNELS)
reason = "Skipped because pythran translation does not work on Mac currently"
nmax_kernel_combinations = [couple if couple[1] != "pythran" else
                            param(couple[0], couple[1],
                                  marks=mark.skipif(sys.platform == "darwin",
                                  reason=reason))
                            for couple in cartesian_product]


@mark.parametrize('nmax, kernel', nmax_kernel_combinations)
def test_nmax_kernel(nmax, kernel):
    """Test a combination of nmax values and kernel versions and compare to
    expected numerical error"""
    tmax = 1.0
    err_ref = {100: 0.04496958454369648,
               1000: 0.018877105158683641}
    error = burgers(tmax, nmax, kernel=kernel)
    assert error == approx(err_ref[nmax])
