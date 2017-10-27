#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that error of the Godunov solver with various kernel versions
"""

from burgers import main, KERNELS
from pytest import mark, approx
import sys

if sys.platform == "darwin":
    # pythran translation does not work on Mac currently
    print("Warning: pythran kernel not tested")
    KERNELS.remove("pythran")


@mark.parametrize('nmax', [100, 1000])
def test_nmax(nmax):
    tmax = 1.0
    err_ref = {100: 0.04496958454369648,
               1000: 0.018877105158683641}
    error = main(tmax, nmax, kernel='python')
    assert error == approx(err_ref[nmax])


@mark.parametrize('kernel', KERNELS)
def test_kernels(kernel):
    tmax = 1.0
    nmax = 100
    error = main(tmax, nmax, kernel=kernel)
    assert error == approx(0.04496958454369648)
