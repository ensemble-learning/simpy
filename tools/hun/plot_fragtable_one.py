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

# plt the reactants


head_label = ['Timestep', r'NH$_3$OH', r'C$_2$O$_2$N$_8$', 'H', r'NH$_{3-x}$OH$_x$', 'H2ON', 'H2N', r'N$_2$O', 'OH', 'HON', 'H$_2$O', 'OC2N6', 'O', 'NH$_3$', 'N$_2$', 'HO2C2N8', 'H3N2', 'H2O2C2N7', 'HOC2N6', 'HO3C2N8', 'H2OC2N6', 'H3O2C2N7', 'H4N2', 'H2O2C2N9', 'H3OC2N6', 'H2O2C2N8', 'HOC2N5', 'H4OC2N6', 'ON3', 'H3O3C2N8', 'H2O3C2N8', 'H4O2C2N7', 'H7O3C2N8', 'HO2C2N7', 'H4N']

Reac = ["H4ON", "O2C2N8"]
fig = plt.figure()
ax = fig.add_subplot(2,2,1)

for i in Reac:
    n = head.index(i)
    ax.plot(data[0]/200.0, data[n], label=head_label[n], lw=3)
ax.set_ylim([0,20])
ax.legend(loc=1)
ax.set_ylabel("Number of Fragments (N)")
ax.set_xlim([13.0, 15])
    
Pro = ["H", "H3ON","ON2", "H2O", "HO", "H3N", "N2"]
ax = fig.add_subplot(2,2,3)

for i in Pro:
    n = head.index(i)
    ax.plot(data[0]/200.0, data[n], label=head_label[n], lw=3)
ax.legend(loc=1)
ax.set_xlabel("Simulation Time (ps)",size="x-large")
ax.set_ylabel("Number of Fragments (N)",size="x-large")
ax.set_xlim([13.0, 15])

data = np.loadtxt("temperature")
ax = fig.add_subplot(2,2,2)
ax.plot(np.arange(len(data))/1000.0,data)
ax.set_xlabel("Simulation Time (ps)", size="x-large")
ax.set_ylabel("Temperature (K)", size="x-large")
ax.set_xlim([13.0, 15])

plt.show()
