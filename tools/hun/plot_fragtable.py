import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import itertools


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
print head, len(head)

# plt the data
cols=plt.cm.Spectral( np.linspace( 0,1,ndata))

cols=itertools.cycle(cols)

fig = plt.figure()
jet = plt.get_cmap('lines')
#cNorm  = colors.Normalize(vmin=0, vmax=(len(head) - 1))
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
nfig = int(math.ceil(ndata/7))

for i in range(1,ndata+1):
    #colorVal = scalarMap.to_rgba(i)
    ax = fig.add_subplot(7,7, i)
    #ax.plot(data[0]/100.0, data[i], label=head[i], lw=4, color=cols.next())
    ax.plot(data[0]/200.0, data[i], label=head[i], lw=3)
    ax.set_title(head[i])
    ax.set_xlim([13.0, np.max(data[0]/200.0)])
    if np.max(data[i]) < 4:
        #ax.plot([9.93, 9.93], [0, 4], ls="--", lw=2)
        #ax.plot([13.76, 13.76], [0, 4], ls="--", lw=2)
        ax.set_ylim([0,4])
    else:
        #ax.plot([9.93, 9.93], [0, np.max(data[i])], ls="--", lw=2)
        #ax.plot([13.76, 13.76], [0, np.max(data[i])], ls="--", lw=2)
        pass
    
plt.show()
