"""
Calculate the diffusion coefficient from msd.
@log:
20141120: add time step
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import poly1d, polyfit

##########################
# 
# msd in A^2 and t in fs (1e-15s)
#
##########################

print "python msd.py [time step (in fs)]"

ts = 1.0

if len(sys.argv) > 1:
    ts = float(sys.argv[1])

toM = 1e-5
toCM = 1e-1
toCM5 = 1e4

data = np.loadtxt("msd.dat")
data = data.transpose()

start = int(0.2*len(data[0]))
end = int(0.8*len(data[0]))

x = data[0] * ts
y = data[1]

fx = x[start:end]
fy = y[start:end]

fit = polyfit(fx, fy, 1)
fit_fn = poly1d(fit)

print "%.4e"%(fit[0]/6*toM)
print "%.4e"%(fit[0]/6*toCM)
print "%.4e 10^-5cm^2/s"%(fit[0]/6*toCM5)

plt.plot(x, y)
plt.plot(fx, fit_fn(fx), '--k', lw=2)
plt.xlabel("Simulation Time (fs)")
plt.ylabel("MSD ($\AA^2$)")
plt.show()

