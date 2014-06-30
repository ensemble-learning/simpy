import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import itertools

timestep = 5 # fs
tops = 1000/timestep

f = open("fragtable", "r")

counter = 0

data = []
head = []

for i in f:
    if counter == 0:
        head = i.strip().split()
    else:
        tokens = i.strip().split()
        if len(tokens) > 1:
            data.append([float(j) for j in tokens])
    counter += 1     

f.close()

data = np.array(data)
data = data.transpose()
ndata = len(head) - 1

nsubx = int(math.sqrt(ndata))
nsuby = int(ndata/nsubx) + 1

# plt the data
cols=plt.cm.Spectral( np.linspace( 0,1,ndata))

cols=itertools.cycle(cols)

fig = plt.figure()
jet = plt.get_cmap('lines')
#cNorm  = colors.Normalize(vmin=0, vmax=(len(head) - 1))
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

for i in range(1,ndata+1):
    #colorVal = scalarMap.to_rgba(i)
    ax = fig.add_subplot(nsubx,nsuby, i)
    #ax.plot(data[0]/100.0, data[i], label=head[i], lw=4, color=cols.next())
    ax.plot(data[0]/tops, data[i], label=head[i], lw=3)
    ax.set_title(head[i], color="red", size="small", weight="bold")
    #ax.set_xlim([13.0, np.max(data[0]/200.0)])
    if np.max(data[i]) < 4:
        #ax.plot([9.93, 9.93], [0, 4], ls="--", lw=2)
        #ax.plot([13.76, 13.76], [0, 4], ls="--", lw=2)
        ax.set_ylim([0,4])
    else:
        #ax.plot([9.93, 9.93], [0, np.max(data[i])], ls="--", lw=2)
        #ax.plot([13.76, 13.76], [0, np.max(data[i])], ls="--", lw=2)
        pass
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
        # specify integer or one of preset strings, e.g.
        #tick.label.set_fontsize('x-small') 
        #tick.label.set_rotation('vertical')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
    # output data
    if not os.path.exists("frag_details"):
        os.mkdir("frag_details")
    os.chdir("frag_details")
    tmp = np.array([data[0], data[i]])
    np.savetxt(head[i] + ".dat", tmp)
    os.chdir("..")

#plot the temperature
if os.path.exists("temperature"):
    tempt = np.loadtxt("temperature")
    ax = fig.add_subplot(nsubx,nsuby, ndata+1)
    ax.plot(data[0]/tops, tempt[1].transpose(), lw=3)
    ax.set_title("Temperature", color="red", size="small", weight="bold")
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
        # specify integer or one of preset strings, e.g.
        #tick.label.set_fontsize('x-small') 
        #tick.label.set_rotation('vertical')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
    
#plot the potential
if os.path.exists("potential"):
    pot = np.loadtxt("potential")
    ax = fig.add_subplot(nsubx,nsuby, ndata+1)
    ax.plot(np.linspace(0,80, len(pot)), pot, lw=3)
    ax.set_title("Potential", color="red", size="small", weight="bold")
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
        # specify integer or one of preset strings, e.g.
        #tick.label.set_fontsize('x-small') 
        #tick.label.set_rotation('vertical')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6) 

fig.subplots_adjust(wspace=0.4, hspace=0.8)

plt.show()

