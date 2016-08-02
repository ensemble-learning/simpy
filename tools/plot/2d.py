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
CS = plt.contour(xi, yi, zi, [-0.80, -0.7, -0.6, -0.5], colors='k')
plt.clabel(CS, inline=1, fontsize=10)

plt.xlim([0.9,2.6])
plt.ylim([1.6,3.6])
plt.xticks([1.0, 1.5, 2.0, 2.5], fontsize=14 )
plt.yticks(fontsize=14 )
plt.xlabel("$\sqrt{r_{O-H_1}^2+r_{O-H_2}^2}$ ($\AA$)", 
           fontsize=24,)
plt.ylabel("$r_{Pt-O}$ ($\AA$)", fontsize=24,)
#plt.scatter(x, y, c=z)
plt.tight_layout()
plt.savefig("fe_2d.eps")
plt.savefig("fe_2d.png", dpi=450)
plt.show()
