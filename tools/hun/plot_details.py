import math
import numpy as np
import matplotlib.pyplot as plt

timestep = 0.10
tops = 1000/timestep


fnames = []
f = open("flist", "r")

for i in f:
    if len(i.strip()) > 0:
        fnames.append(i.strip().split())

f.close()


nx = int(math.sqrt(len(fnames)))
ny = int(len(fnames)/nx) + 1

fig = plt.figure()

for i in range(len(fnames)):
    fname = fnames[i][0]
    molname = fnames[i][1]
    if fname == "blank":
        pass
    else:
        ax = fig.add_subplot(nx, ny, i+1)
        data = np.loadtxt(fname)
        ax.plot(data[0]/tops, data[1], lw=2)
        ax.set_title(molname, color="red", size="small", weight="bold")
        ymax = 5*(int(np.max(data[1])/5) + 1)
        ax.set_ylim([0,ymax])

fig.subplots_adjust(wspace=0.4, hspace=0.8)
plt.show()

