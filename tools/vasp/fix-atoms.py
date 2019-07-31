import shutil, sys

n0 = 9
fix_atoms = []

g1 = range(29)
for i in g1:
    fix_atoms.append(i)

"""
g2 = [2, 10, 11, 19]
for i in g2:
    fix_atoms.append(i)
"""

f = open("POSCAR", "r")
poscar = f.readlines()
f.close()

# the first line poscar[9]
shutil.copy("POSCAR", "POSCAR.0")

for i in fix_atoms:
    nid = n0 + i
    poscar[nid] = poscar[nid].replace("T", "F")

o = open("POSCAR_new", "w")
for i in poscar:
    o.write(i)
o.close()

