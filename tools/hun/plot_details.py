import math
import numpy as np
import matplotlib.pyplot as plt
import os

timestep = 0.25
tops = 1/timestep

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

if os.path.exists("temperature"):
    nplot += 1

if os.path.exists("potential"):
    nplot += 1

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
        ax.plot(data[0]/4000, data[1], lw=2, label=molname)
        #ax.set_title(molname, color="red", size="small", weight="bold")
        #ymax = 5*(int(np.max(data[1])/5) + 1)
        #ax.set_ylim([0,ymax])

ax.legend(loc=1)
ax.set_xlabel("Simulation Time (ps)")
#fig.subplots_adjust(wspace=0.4, hspace=0.8)
#plot the temperature if any
if os.path.exists("temperature"):
    data = np.loadtxt("temperature")
    ax = fig.add_subplot(nplot, 1, n + 1)
    ax.plot(data)
    ax.set_xlabel("Temperature (K)")

#plot the potential energy if any
if os.path.exists("potential"):
    data = np.loadtxt("potential")
    ax = fig.add_subplot(nplot, 1, n + 2)
    ax.plot(data)
    ax.set_xlabel("Potential Energy (ev/mol)")

plt.show()

