import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("all.dat")
data = data.transpose()

hist_data = []
n_total = 0
for i in range(len(data[0])):
    ndata = int(data[1][i])
    n_total += ndata
    for j in range(ndata):
        hist_data.append(data[0][i])

hist, bin_edges = np.histogram(hist_data, 12)
dx = bin_edges[1] - bin_edges[0]
o = open("all-hist.dat", "w")
for i in range(len(hist)):
    o.write("%12.8f%12.4f\n"%(bin_edges[i]+0.5*dx, 
                              hist[i]*36.0/n_total))
o.close()

plt.plot(bin_edges[:-1], hist*36.0/n_total)
plt.show()


    
