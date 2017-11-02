#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test if results are identical to the reference data from C version
"""

from stvenant import wL, wR, xi, load_file
from stvenant import riemann_loop, KERNELS, C_program
from pytest import fixture, mark, param
from numpy.testing import assert_allclose, assert_array_equal
import sys


@fixture
def ref_data():
    """A fixture that provides the reference data"""
    return load_file("plotriem_nx1000.dat")


def test_xi(ref_data):
    """Test if space vector are identical"""
    xi_ref, h_ref, u_ref = ref_data
    assert_allclose(xi, xi_ref, atol=1e-15)


@mark.parametrize('kernel', KERNELS)
def test_loop(ref_data, kernel):
    """Test xi-loop riemann_loop for various implementations of the
    riemmann solver function"""
    xi_ref, h_ref, u_ref = ref_data
    w = riemann_loop(kernel, wL, wR, xi)
    assert_allclose(w[:, 0], h_ref, atol=1e-6)
    assert_allclose(w[:, 1]/w[:, 0], u_ref, atol=1e-6)


def test_C(ref_data):
    """Test if stvenant.exe produces the same data output"""
    x_ref, h_ref, u_ref = ref_data
    C_program("./stvenant.exe")
    x, h, u = load_file("plotriem")
    assert_array_equal(x, x_ref)
    assert_array_equal(h, h_ref)
    assert_array_equal(u, u_ref)
