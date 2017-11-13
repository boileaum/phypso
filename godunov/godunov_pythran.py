#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Godunov solver
"""

from godunov_python import Burgers as Burgers_python


class Burgers(Burgers_python):

    def __init__(self):
        from riemann_pythran import riemann_burgers
        self.riemman = riemann_burgers
