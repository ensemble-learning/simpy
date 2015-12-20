"""
Plot the neb results
"""
import numpy as np
import matplotlib.pyplot as plt

neb = []
f = open("log.lammps", "r")

for i in f:
    if i.strip().startswith("Step"):
        break

for i in f:
    if i.strip().startswith("Climbing"):
        break
    else:
        tokens = i.strip().split()
        data = tokens[9:]
        tmp = [[],[]]
        for j in range(0, len(data), 2):
            lamda = float(data[j])
            de = float(data[j+1])
            tmp[0].append(lamda)
            tmp[1].append(de)
        neb.append(tmp)
    
for i in f:
    if i.strip().startswith("Step"):
        break

for i in f:
    tokens = i.strip().split()
    data = tokens[9:]
    tmp = [[],[]]
    for j in range(0, len(data), 2):
        lamda = float(data[j])
        de = float(data[j+1])
        tmp[0].append(lamda)
        if j == 0:
            e0 = de
        tmp[1].append(de - e0)
    neb.append(tmp)

plt.plot(neb[-1][0], neb[-1][1], "-o")
plt.show()
