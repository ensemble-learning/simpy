#!/usr/bin/env python

"""
convert the exchange current density to reaction barrier.
ref: acscatal.7b03142
note: tau = 5e14 for RuO2
note: tau = 2.6e15 for Pt(111)
"""

import sys, math
#constants
kb = 1.38064852e-23 # iBoltzmann constant (J K-1)
e = 1.6021766208e-19   # Elementary charge (C)
h = 6.62607004e-34  # Planck's constant (m2 kg s-1)
na = 6.022140857e23 # Avogadro constant (mol-1)
ev2jmol = 96485.000 # ev to jmol-1

#variables
T = 300             # Temperautre (K)
z = 1               # number of electrons transferred in the overall reaction (n)
tau = 2.6e15          # number of electrocatalyst's active sites per surface area (cm-2)
j0 = 1e-6           # exchange current density (A cm-2)
G = 0.89            # TS free energy (eV)

def test():
    a1 = kb*T*z*e*tau/h
    a2 = math.exp(-G*ev2jmol/(kb*na*T))
    j0_log = math.log10(a1*a2)
    print(a1, a2, j0_log)

if 0:
    test()

def j2g():
    if len(sys.argv) > 1:
        j0_log = float(sys.argv[1])
        j0 = math.pow(10, j0_log)
        a1 = kb*T*z*e*tau/h
        a2 = j0/a1
        a3 = math.log(a2)
        G = -a3*kb*na*T/ev2jmol
        print("At %.2f K, with active sites density of %.2e cm-2"%(T, tau))
        print("A current density of %.2f (log10 Acm-2) leads to:"%j0_log)
        print("dG# = %.2f"%G, "eV")
    else:
        print("j2g log10(j0) [sites density (cm-1)] [Temperature]")

j2g()
