import sys
import numpy as np
import matplotlib.pyplot as plt

infile = "data.dat"
if len(sys.argv) > 1:
    infile = sys.argv[1]

data = np.loadtxt(infile)
hist, edge = np.histogram(data, bins=51)
dx = edge[1] - edge[0]
x = edge[:-1] + 0.5*dx

plt.bar(x, hist, width=dx*0.8) 
#plt.show()
plt.savefig("hist.png")

o = open("hist.dat", "w")
for i in range(len(x)):
    o.write("%12.4f%12.4f\n"%(x[i], hist[i]))
o.close()
