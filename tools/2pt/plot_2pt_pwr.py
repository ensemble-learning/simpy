import numpy as np
import matplotlib.pyplot as plt

data = []

f = open("output.pwr", "r")

counter = 0

for i in f:
    tokens = i.strip().split()
    if counter < 0:
        pass
    elif counter == 2:
        title = tokens
        n = len(title)
    elif counter > 2:
        if len(tokens) == n:
            data.append([float(j) for j in tokens])
    counter += 1

f.close()

data = np.array(data)
data = data.transpose()

for i in range(1, len(data)):
    if i in [5,6,7,8,9]:
        plt.plot(data[0], data[i], label=title[i], lw=2)

plt.xlabel("Frequency (cm$^{-1}$)", size="x-large")
plt.ylabel("N", size="x-large")
plt.legend()
plt.show()
