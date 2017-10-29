#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test if results are identical to reference data from C version
"""

from stvenant import wL, wR, xi, load_file, riemann_loop, riemann_loop_C
from pytest import fixture


@fixture
def ref_data():
    """A fixture that provides the reference data"""
    return load_file("plotriem_nx1000.dat")


def test_xi(ref_data):
    """Test if space vector are identical"""
    x_ref, h_ref, u_ref = ref_data
    assert xi.all() == x_ref.all()


def test_loop_python(ref_data):
    """Test full python version"""
    x_ref, h_ref, u_ref = ref_data
    w = riemann_loop(wL, wR, xi)
    assert w[:, 0].all() == h_ref.all()
    assert w[:, 1].all() == u_ref.all()


def test_loop_C(ref_data):
    """Test python calling a C kernel"""
    x_ref, h_ref, u_ref = ref_data
    w = riemann_loop_C(wL, wR, xi)
    assert w[:, 0].all() == h_ref.all()
    assert w[:, 1].all() == u_ref.all()
