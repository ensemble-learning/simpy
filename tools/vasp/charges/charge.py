"""
Parse the charges obtained from bader charge analysis.
"""

import numpy as np
import os

def read_charges():
    """
    Read bader charges from ACF.dat
    """

    # Read bader charges
    lines = []

    f = open("ACF.dat", "r")

    for i in f:
        if i.strip().startswith("--"):
            break

    for i in f:
        if i.strip().startswith("--"):
            break
        else:
            lines.append(i.strip().split())
    f.close()
    return lines

def read_ndx():
    """ Read index
    """
    # Reading index
    ndx = []
    f = open("index.ndx", "r")

    for i in f:
        tokens = i.strip().split()
        for j in tokens:
            ndx.append(int(j))

    f.close()
    return ndx

def read_nele():
    """ Read nelectron
    """
    # Reading nelectron
    nele = []
    f = open("atoms.dat", "r")

    for i in f:
        tokens = i.strip().split()
        for j in tokens:
            nele.append(int(j))

    f.close()
    return nele


def sum_charges(lines, ndx, nele):
    """ Calculate the charges
    """

    charges = []

    for i in range(len(ndx)):
        tokens = lines[ndx[i]-1]
        #print tokens
        q1 = float(tokens[4])
        q2 = nele[i] - q1
        #print q1, q2
        charges.append(q2)

    charges = np.array(charges)

    print np.sum(charges)

def main():
    
    lines = read_charges()
    ndx = read_ndx()
    nele = read_nele()
    sum_charges(lines, ndx, nele)

if __name__ == "__main__":
    main()

