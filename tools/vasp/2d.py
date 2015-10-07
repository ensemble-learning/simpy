import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

# Generate data:
data = np.loadtxt("fes.dat")
data = data.transpose()

x = data[0]
y = data[1]
z = data[2]

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
                      extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.xlabel("CV 1", size="x-large")
plt.ylabel("CV 2", size="x-large")
CS = plt.contour(xi, yi, zi, 6, colors='k')
plt.clabel(CS, inline=1, fontsize=10)
#plt.scatter(x, y, c=z)
plt.savefig("fe_2d.eps")
plt.savefig("fe_2d.png", dpi=450)
plt.show()
