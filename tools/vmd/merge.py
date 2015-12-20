import numpy as np

r1 = np.loadtxt("r1.dat")
r1 = r1.transpose()

r2 = np.loadtxt("r2.dat")
r2 = r2.transpose()

r = []

for i in range(len(r1[1])):
    if r1[1][i] < r2[1][i]:
        val = r1[1][i]
    else:
        val = r2[1][i]
    r.append(val)

o = open("r.dat", "w")
for i in range(len(r)):
    o.write("%20.8f\n"%r[i])
o.close()

