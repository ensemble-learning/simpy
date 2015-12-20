import sys
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

x0 = -5.05
x0 = -0.114047310
x1 = -0.014047310

x0 = -0.554047310
x1 = -0.014047310

print (x1 - x0) * 27.0

fp = sys.argv[1]
data = np.loadtxt(fp, skiprows=1)
data = data.transpose()

x, y = [], []
for i in range(len(data[0])):
    if data[0][i] <= x1 and data[0][i] > x0:
        x.append(data[0][i])
        y.append(data[1][i])
#print len(x), len(y)

y_int = integrate.cumtrapz(y, x, initial=0)

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(x, y)

ax = fig.add_subplot(2,1,2)
ax.plot(x, y_int)
print y_int[-1]

plt.show()
