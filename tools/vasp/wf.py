#!/usr/bin/env python

"""
Calculate work function from VASP output.
"""

import os

def get_efermi():
    f = open("OUTCAR", "r")
    for i in f:
        if i.strip().startswith("E-fermi"):
            tokens = i.strip().split()
            efermi = float(tokens[2])
    return efermi

def get_eq():
    n0 = 120
    n1 = 130
    f = open("vplanar.txt", "r")
    n = 0
    eq = 0.0
    for i in f:
        if n >n0 and n <= n1:
            tokens = i.strip().split()
            eq += float(tokens[1])
        n += 1
    eq = eq/(n1-n0)
    return eq
        
        
assert os.path.exists("OUTCAR")
assert os.path.exists("vplanar.txt")
efermi = get_efermi()
eq = get_eq()
wf = eq - efermi
print "%.2f"%wf

