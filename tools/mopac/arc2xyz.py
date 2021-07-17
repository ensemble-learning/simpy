#! /usr/bin/env python
"""
Get the geo from the mopac (.arc)
"""

import sys
import math

infile = sys.argv[1]

f = open(infile, "r")
for i in f:
    if i.strip().startswith("FINAL GEOMETRY OBTAINED"):
        break

for i in f:
    tokens = i.strip()
    if len(tokens) == 0:
        break

lines = []
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 7:
        lines.append(tokens)

f.close()

o = open("geo_end.xyz", "w")
o.write("%d\n\n"%len(lines))
for i in lines:
    print(i)
    ele = i[0]
    x = float(i[1]) 
    y = float(i[3]) 
    z = float(i[5])
    o.write("%s\t%.6f\t%.6f\t%.6f\n"%(ele, x, y, z))
o.close()
