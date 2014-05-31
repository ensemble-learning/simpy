import math
import numpy as np
import matplotlib.pyplot as plt

timestep = 0.10
tops = 100/timestep

nplot = 1

fnames = []
f = open("flist", "r")

for i in f:
    if len(i.strip()) > 0:
        if i.startswith("#"):
            pass
        else:
            if "blank" in i:
                nplot += 1
            fnames.append(i.strip().split())

f.close()

fig = plt.figure()

n = 1
for i in range(len(fnames)):
    fname = fnames[i][0]
    molname = fnames[i][1]
    if fname == "blank":
        if n > 1:
            ax.legend(loc=2)
        else:
            ax.legend()
        n += 1
    else:
        ax = fig.add_subplot(nplot, 1, n)
        data = np.loadtxt(fname)
        ax.plot(data[0]/tops, data[1], lw=2, label=molname)
        #ax.set_title(molname, color="red", size="small", weight="bold")
        #ymax = 5*(int(np.max(data[1])/5) + 1)
        #ax.set_ylim([0,ymax])

ax.legend(loc=2)
ax.set_xlabel("Simulation Time (ps)")
#fig.subplots_adjust(wspace=0.4, hspace=0.8)
plt.show()

