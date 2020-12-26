import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

fname = 'xrd.dat'

data = np.loadtxt(fname, skiprows=4)
data = data.transpose()

fig = plt.figure(figsize=(4, 3), dpi=300)
ax = fig.add_subplot(111)

ax.plot(data[1], data[2], color='black', lw=0.5)
#ax.scatter(data[1], data[2], 8, color='black')

ax.set_xlim([5.0, 50.0])
ax.set_xlabel(r'2$\theta$ (deg)', fontsize=14)
ax.set_ylabel('Intensity', fontsize=14)

plt.tight_layout()
plt.savefig('xrd.svg')
