"""Plot the pressure, potential energy and temperature of the simulation
"""

import matplotlib.pyplot as plt
import numpy as np

step = []
press = []
pot = []
tempt = []

f = open("md.out", "r")

for i in f:
    tokens = i.strip().split(":")
    if "MD step" in i:
        val = tokens[1].split()[0]
        step.append(int(val))
    if "Pressure" in i:
        val = tokens[1].split()[2]
        press.append(float(val))
    if "Potential Energy" in i:
        val = tokens[1].split()[2]
        pot.append(float(val))
    if "Temperature" in i:
        val = tokens[1].split()[2]
        tempt.append(float(val))
        
f.close()

step = np.array(step)
press = np.array(press)
pot = np.array(pot)
tempt = np.array(tempt)

fig = plt.figure()

ax = fig.add_subplot(3,1,1)
ax.set_title("Pressure")
ax.set_ylabel("Pa")
ax.plot(step, press)
np.savetxt("pressure", np.array([step, press]))

ax = fig.add_subplot(3,1,2)
ax.set_title("Potential Energy")
ax.set_ylabel("ev")
ax.plot(step, pot)
np.savetxt("potential", np.array([step, pot]))

ax = fig.add_subplot(3,1,3)
ax.set_title("Temperature")
ax.set_ylabel("K")
ax.plot(step, tempt)
np.savetxt("temperature", np.array([step, tempt]))

plt.show()

