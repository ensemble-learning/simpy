import os
from fort import Fort73
import numpy as np
import matplotlib.pyplot as plt

os.chdir("/home/tao/Documents/wag/training/ca")
os.system("./exe.sh")
def addPlot(filename, tag):
    data = Fort73(filename)
    #ebond = np.array(data.ebond) - data.ebond[11]
    #eatom = np.array(data.eatom) - data.eatom[11]
    #evdw = np.array(data.evdw) -data.evdw[11]

    ebond = np.array(data.ebond)
    eatom = np.array(data.eatom)
    evdw = np.array(data.evdw)

    #plt.plot(ebond, label=tag+' ebond')
    #plt.plot(eatom, label=tag+' eatom')
    #plt.plot(evdw, label=tag+' evdw')
    plt.plot(ebond+eatom+evdw, label=tag)

addPlot('in/fort.73', tag='increase')
addPlot('de/fort.73', tag='decrease')
addPlot('back.fort.73', tag='old')

qm = [ -0.946, -0.768, -0.611, -0.474, -0.350, -0.241,\
-0.142, -0.077, -0.034, -0.009, -0.001, -0.011,\
-0.048, -0.112, -0.200, -0.315, -0.437, -0.598,\
-0.791, -1.012, -1.255,]
qm = np.array(qm)
qm = -72.48 - qm
plt.plot(qm, 'o', label="qm")
plt.legend()
plt.show()
