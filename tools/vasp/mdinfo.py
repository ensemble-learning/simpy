"""
get MD information from VASP QMD calculation
"""

import sys
import numpy as np

if len(sys.argv) < 2:
    t0 = 0
else:
    t0 = int(sys.argv[1])

steps = []
pote = []
potf = []
T = []

fname = "OSZICAR"
f = open(fname, "r")

for i in f:
    if "=" in i:
        tokens = i.strip().split()
        steps.append(int(tokens[0]))
        T.append(float(tokens[2]))
        potf.append(float(tokens[6]))

f.close()

o = open("potential", "w")
for i in potf:
    o.write("%.4f\n"%i)
o.close()

o = open("temperature", "w")
for i in T:
    o.write("%.4f\n"%i)
o.close()
    
flag_plot = 0
if flag_plot:
    import matplotlib.pyplot as plt

    plt.plot(pote)
    plt.show()

potf = np.array(potf)
print(np.average(potf[t0:]))

