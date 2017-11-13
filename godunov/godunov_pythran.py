#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver
"""

from godunov_python import Burgers as Burgers_python


class Burgers(Burgers_python):

    def __init__(self, nmax, tmax):
        super().__init__(nmax, tmax)
        from riemann_pythran import numflux
        self.numflux = numflux
