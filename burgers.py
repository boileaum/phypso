#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 18:48:24 2017

@author: boileau
"""

import argparse
import godunov
import matplotlib.pyplot as plt


def main(tmax):

    xm, wn = godunov.timeloop(tmax)
    wex = godunov.sol_exact(tmax)
    with open('exact.dat', 'w') as f_exact:
        with open('godu.dat', 'w') as f_godu:
            for i in range(1, godunov.NMAX+1):
                f_exact.write("{} {}\n".format(xm[i], wex[i]))
                f_godu.write("{} {}\n".format(xm[i], wn[i]))

    plt.plot(xm[1:-1], wex[1:-1], 'r+', label='exact')
    plt.plot(xm[1:-1], wn[1:-1], 'k-', label='godunov')
    plt.legend()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Solve Burgers problem")
    parser.add_argument(dest='tmax', nargs=1, metavar="final_time",
                        help="Simulation final time")
    args = parser.parse_args()
    tmax = float(args.tmax[0])
    print(tmax)
    main(tmax)
