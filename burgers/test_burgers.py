#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that numerical error of the Godunov solver with various kernel versions
"""

from burgers import burgers, KERNELS
from pytest import mark, approx
import itertools
import sys

if sys.platform == "darwin":
    # pythran translation does not work on Mac currently
    print("Warning: pythran kernel not tested")
    KERNELS.remove("pythran")


@mark.parametrize('nmax, kernel', itertools.product([100, 1000], KERNELS))
def test_nmax_kernel(nmax, kernel):
    """Test a combination of nmax values and kernel versions and compare to
    expected numerical error"""
    tmax = 1.0
    err_ref = {100: 0.04496958454369648,
               1000: 0.018877105158683641}
    error = burgers(tmax, nmax, kernel=kernel)
    assert error == approx(err_ref[nmax])
