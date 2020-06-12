#!/usr/bin/env python

"""
usage: k2c.py tafel_slope
convert the reaction barriers to current density
ref: acscatal.7b03142
"""

import sys, math
import numpy as np

#constants
kb = 1.38064852e-23 # iBoltzmann constant (J K-1)
e = 1.6021766208e-19   # Elementary charge (C)
h = 6.62607004e-34  # Planck's constant (m2 kg s-1)
na = 6.022140857e23 # Avogadro constant (mol-1)
ev2jmol = 96485.000 # ev to jmol-1

def get_current(eta, T, tau, a1, a2, gact1, gact2, dgtd1):
    c0 = 2*kb*T*e*tau/h
    n1 = ((a1 + a2)*eta*e)/(kb*T)
    n1 = math.exp(n1)

    n2 = (-(2 - a1 - a2)*eta*e)/(kb*T)
    n2 = math.exp(n2)

    d1 = (gact1*ev2jmol/na + a2*eta*e)/(kb*T)
    d1 = math.exp(d1)

    d2 = (gact2*ev2jmol/na - dgtd1*ev2jmol/na + a1*eta*e)/(kb*T)
    d2 = math.exp(d2)

    d3 = (gact1*ev2jmol/na - dgtd1*ev2jmol/na - (1-a2)*eta*e)/(kb*T)
    d3 = math.exp(d3)

    d4 = (gact2*ev2jmol/na - (1-a1)*eta*e)/(kb*T)
    d4 = math.exp(d4)

    ju = c0*(n1-n2)/(d1+d2+d3+d4)
    return ju, math.log10(ju)

def test_cer():
    # variables for CER
    T = 300.0         # Temperautre (K)
    tau = 5e14        # number of electrocatalyst's active sites per surface area (cm-2)
    a1 = 0.69         # symmetry factors
    a2 = 0.64         # symmetry factors
    gact1 = 0.77      # TS free energy
    gact2 = 0.89      # TS free energy
    dgtd1 = 0.3       # free energy difference of reaction intermediate S-X
    eta = np.arange(0.02, 0.21, 0.01)
    for i in eta:
        ju, ju_log10 = get_current(i, T, tau, a1, a2, gact1, gact2, dgtd1)
        print("%8.4f%10.6f%10.4f"%(i, ju, ju_log10))

def test_her():
    # variables for HER
    T = 300.0         # Temperautre (K)
    tau = 2.6e15        # number of electrocatalyst's active sites per surface area (cm-2)
    a1 = 0.405        # symmetry factors
    a2 = 0.74         # symmetry factors
    gact1 = 0.665     # TS free energy
    gact2 = 0.75      # TS free energy
    dgtd1 = 0.35      # free energy difference of reaction intermediate S-X
    eta = np.arange(0.02, 0.21, 0.01)
    for i in eta:
        ju, ju_log10 = get_current(i, T, tau, a1, a2, gact1, gact2, dgtd1)
        print("%8.4f%10.6f%10.4f"%(i, ju, ju_log10))

def main():
    eta = np.arange(0.01, 0.12, 0.01)
    for i in eta:
        ju, ju_log10 = get_current(i, T, tau, a1, a2, gact1, gact2, dgtd1)
        print("%8.3f%10.6f%10.4f"%(i, ju, ju_log10))

if __name__ == "__main__":
    #test_cer()
    test_her()
    #main()
