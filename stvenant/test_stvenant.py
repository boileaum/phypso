#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test if results are identical to the reference data from C version
"""

from stvenant import wL, wR, xi, load_file
from stvenant import riemann_loop, KERNELS
from pytest import fixture, mark, param
import sys


@fixture
def ref_data():
    """A fixture that provides the reference data"""
    return load_file("plotriem_nx1000.dat")


def test_xi(ref_data):
    """Test if space vector are identical"""
    x_ref, h_ref, u_ref = ref_data
    assert xi.all() == x_ref.all()


reason = "Skipped because pythran translation does not work on Mac currently"
kernels = [kernel if kernel != "pythran"
           else param(kernel, marks=mark.skipif(sys.platform == "darwin",
                      reason=reason)) for kernel in KERNELS]


@mark.parametrize('kernel', kernels)
def test_loop(ref_data, kernel):
    """Test xi-loop riemann_loop for various implementations of the
    riemmann solver function"""
    x_ref, h_ref, u_ref = ref_data
    w = riemann_loop(wL, wR, xi, kernel)
    assert w[:, 0].all() == h_ref.all()
    assert w[:, 1].all() == u_ref.all()
