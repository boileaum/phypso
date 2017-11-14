#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver
"""

from godunov_python import Burgers as Burgers_python
from godunov_python import StVenant as StVenant_python


class Burgers(Burgers_python):

    def __init__(self, nmax, tmax):
        super().__init__(nmax, tmax)
        from riemann_burgers_pythran import numflux
        self.numflux = numflux


class StVenant(StVenant_python):

    def __init__(self, nmax, tmax):
        super().__init__(nmax, tmax)
        from riemann_stvenant_pythran import riemann as \
            riemann_stvenant_pythran
        self.riemann = riemann_stvenant_pythran
