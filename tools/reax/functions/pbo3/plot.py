import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data1")
data = data.transpose()

plt.plot(data[0], data[1], lw=2, color="blue")
plt.plot(data[0], data[2], lw=2, color="blue", ls="--")

data = np.loadtxt("data2")
data = data.transpose()

plt.plot(data[0], data[1], lw=2, color="green")
plt.plot(data[0], data[2], lw=2, color="green", ls="--")

plt.show()
