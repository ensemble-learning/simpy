import sys
import numpy as np
from scipy import integrate

#datafile = sys.argv[1]
datafile = "fe.dat"

data = np.loadtxt(datafile)
data = data.transpose()

x = data[0]
y = data[1]

y_int = integrate.cumtrapz(y, x, initial=0)

o = open("fe-int.dat", "w")
for i in range(len(x)):
    o.write("%12.4f%12.4f\n"%(x[i], y[i]))

