#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test if results are identical to the reference data from C version
"""

from stvenant import wL, wR, xi, load_file
from stvenant import riemann_loop, riemann_C, riemann_python, riemann_cython
from pytest import fixture, mark


@fixture
def ref_data():
    """A fixture that provides the reference data"""
    return load_file("plotriem_nx1000.dat")


def test_xi(ref_data):
    """Test if space vector are identical"""
    x_ref, h_ref, u_ref = ref_data
    assert xi.all() == x_ref.all()


@mark.parametrize('riemann_func', [riemann_python, riemann_C, riemann_cython])
def test_loop(ref_data, riemann_func):
    """Test xi-loop riemann_loop for various implementations of the
    riemmann solver function"""
    x_ref, h_ref, u_ref = ref_data
    w = riemann_loop(riemann_func, wL, wR, xi)
    assert w[:, 0].all() == h_ref.all()
    assert w[:, 1].all() == u_ref.all()
