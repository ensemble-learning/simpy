#!/usr/bin/env python

import os, math, re
from periodictable import formula
pdbfile = "water-0056.pdb"
os.system("python2 ~/soft/simpy/lib/e_2_pdb.py %s -c"%pdbfile)

amc = 1.66053904e-27
a2m = 1e30

formula_box = {}
pattern = re.compile(r'(\D+)(\d*)')

f = open(pdbfile, "r")
for i in f:
    if i.strip().startswith("CRYST1"):
        a, b, c, alpha, beta, gamma = [float(ii) for ii in i.strip().split()[1:7]]
    if i.strip().startswith("HETATM") or i.strip().startswith("ATOM"):
        tokens = i.strip().split()[2]
        match = pattern.match(tokens)
        if match:
            ele = match.group(1)
            if len(ele) > 1:
                ele = ele[0].upper() + ele[1].lower()
        if ele in formula_box.keys():
            formula_box[ele] += 1
        else:
            formula_box[ele] = 1
formula_box_txt = ''
for i in formula_box.keys():
    formula_box_txt += "%s%d"%(i, formula_box[i])
box = formula(formula_box_txt)

alpha = alpha*math.pi/180.0
beta = beta*math.pi/180.0
gamma = gamma*math.pi/180.0

volume = a*b*c*math.sqrt(1 - math.cos(alpha)*math.cos(alpha) -
                math.cos(beta)*math.cos(beta) -
                math.cos(gamma)*math.cos(gamma) +
                2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))

box.density = box.mass/volume*amc*a2m/1000.0

o = open("UFF.atoms", "w")
o.write("""Pt       3.000
O    2.0
H    2.0
EOF
""")
o.close()

o = open("input.dat", "w")
o.write("""UFF.atoms
out.music
1.400
5000
%-12.4f%-12.4f%-12.4f
%-12.4f


! file containing the atom types + diameters
! file containing the coordinates
! probe size in A
! number of insertions
! length of unitcell
! crystal density in g / cm3
"""%(a, b, c, box.density))
o.close()
