"""
Calculate the diffusion coefficient from msd.
"""
import numpy as np
import matplotlib.pyplot as plt
from pylab import poly1d, polyfit

##########################
# 
# msd in A^2 and t in fs (1e-15s)
#
##########################

toM = 1e-5
toCM = 1e-1
toCM5 = 1e4

data = np.loadtxt("msd.dat")
data = data.transpose()

start = int(0.2*len(data[0]))
end = int(0.8*len(data[0]))

x = data[0][start:end]
y = data[1][start:end]

fit = polyfit(x, y, 1)
fit_fn = poly1d(fit)

print fit[0]/6*toM
print fit[0]/6*toCM
print fit[0]/6*toCM5

plt.plot(data[0], data[1], data[0],fit_fn(data[0]), '--k')
plt.show()

