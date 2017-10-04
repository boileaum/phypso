#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 18:48:24 2017

@author: boileau
"""

import godunov
import matplotlib.pyplot as plt


def main():

    xm, wn = godunov.timeloop()
    wex = godunov.sol_exact()
    with open('exact.dat', 'w') as f_exact:
        with open('godu.dat', 'w') as f_godu:
            for i in range(1, godunov.NMAX+1):
                f_exact.write("{} {}\n".format(xm[i], wex[i]))
                f_godu.write("{} {}\n".format(xm[i], wn[i]))

    plt.plot(xm[1:-1], wex[1:-1], 'r+', label='exact')
    plt.plot(xm[1:-1], wn[1:-1], 'k-', label='godunov')
    plt.legend()

if __name__ == '__main__':
    main()
