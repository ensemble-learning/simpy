#! /usr/bin/env python

"""
Conver the HILLSPOT file (VASP) to HILLS (plumed)
"""

hills = []
f = open("HILLSPOT", "r")
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 3:
        hills.append([float(j) for j in tokens])
f.close()

o = open("HILLS", "w")
o.write("#! FIELDS time d1 sigma_d1 height biasf\n")
o.write("#! SET multivariate false\n")
n = 1
for i in hills:
    o.write("%8d %12.6f %12.6f %12.6f %4d\n"%(n, i[0], i[2], i[1], 1))
    n += 1
o.close()

