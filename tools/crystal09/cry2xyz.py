#! /usr/bin/env python

import sys

infile = sys.argv[1]

f = open(infile, "r")

for i in f:
    if i.strip().startswith("FINAL OPTIMIZED GEOMETRY"):
        break

n = 0

for i in f:
    if n == 4:
        break
    n += 1

atoms = []
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 0:
        break
    else:
        atp = tokens[3]
        x = tokens[4]
        y = tokens[5]
        z = tokens[6]
        if atp == "XX":
            pass
        else:
            atoms.append([atp, x, y, z])
        
f.close()

o = open("output.xyz", "w")
o.write("%d\n"%(len(atoms)))
o.write("\n")
for i in atoms:
    o.write("%s\t%s\t%s\t%s\n"%(i[0], i[1], i[2], i[3]))
o.close()
