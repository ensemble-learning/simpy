"""
Plot the fragments in three subfigures: reactants, productions
and intermediates.
Need a flist file to organize the file names.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

class Frame():
    def __init__(self,):
        self.fnames = []
        self.data = []

# define color map
COLOR = ["#000000", "#5DA5DA", "#FAA43A", "#60BD68", 
         "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", 
         "#F15854", "#4D4D4D"]
nc = 0 # label for color map

# simulation details
timestep = 0.25
tops = 1/timestep

# Read the file names from flist
frames = []
frame = Frame()

f = open("flist", "r")

for i in f:
    if len(i.strip()) > 0:
        if i.startswith("#"):
            pass
        else:
            if "blank" in i:
                frames.append(frame)
                frame = Frame()
            else:
                frame.fnames.append(i.strip().split())
frames.append(frame)
f.close()

# Read the data
nplot = 1
for i in range(len(frames)):
    for j in frames[i].fnames:
        fname = j[0]
        tmp = np.loadtxt(fname, delimiter=",")
        tmp = tmp.transpose()
        frames[i].data.append(tmp)

# determine the space of subplots
if len(frames) > 2:
    nplot = 2 + len(frames[2].data)
    ns = len(frames[2].data)

fig = plt.figure()
gs = gridspec.GridSpec(ns*3, 1)

# plot the first figure
ax = fig.add_subplot(gs[0:ns,0])
for i in range(len(frames[0].data)):
    ax.plot(frames[0].data[i][0], frames[0].data[i][1], label=frames[0].fnames[i][1], color=COLOR[nc])
    nc += 1
ax.xaxis.set_ticklabels([])
ax.set_ylim([0,18])
ax.set_ylabel("N$_{rectants}$")
ax.legend()

# plot the second figure
ax = fig.add_subplot(gs[ns:2*ns,0])
for i in range(len(frames[1].data)):
    ax.plot(frames[1].data[i][0], frames[1].data[i][1], label=frames[1].fnames[i][1], color=COLOR[nc])
    nc += 1
ax.xaxis.set_ticklabels([])
ax.set_ylim([0,17])
ax.set_ylabel("N$_{products}$")
ax.legend(loc=2, ncol=3)

# plot the rest figures
counter = 0
for i in range(len(frames[2].data)):
    ax = fig.add_subplot(gs[2*ns + counter, 0])
    for j in frames[2].data:
        ax.plot(np.linspace(1200, 3000, len(j[0])), j[1], lw=2, color="gray", alpha=0.25)
    ax.plot(np.linspace(1200, 3000, len(frames[2].data[i][0])), frames[2].data[i][1], lw=2, color="black")
    ax.set_frame_on(False)
    ax.get_xaxis().tick_bottom()
    #ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_visible(False)
    molname = frames[2].fnames[i][1]
    #ax.set_ylabel(molname, rotation=0, size=14)
    ax.text(0.05, 0.5, molname, 
            horizontalalignment='center',
             verticalalignment='center',
            transform = ax.transAxes)
    counter += 1
ax.axes.get_yaxis().set_ticks([])
ax.axes.get_xaxis().set_visible(True)
ax.set_xlabel("Temperature (K)")

#fig.subplots_adjust(hspace=0.2)
plt.savefig("rxn.eps")
plt.savefig("rxn.png", dpi=600)
plt.show()

