import numpy as np
import matplotlib.pyplot as plt

data = []

f = open("result.pwr", "r")

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

fig = plt.figure()

ax = fig.add_subplot(3,2,1)
ax.plot(data[0], data[9], label="s($v$)") # total
ax.plot(data[0], data[5], label="s$_{trn}$($v$)") # translation
ax.plot(data[0], data[6], label="s$_{rot}$($v$)") # rotation
ax.set_xlim([0,1000])
ax.set_xlabel("$v$ (cm$^{-1}$)", size="x-large")
ax.set_ylabel("s($v$) (cm)", size="x-large")
ax.legend()

ax = fig.add_subplot(3,2,3)
ax.plot(data[0], data[5], label="s$_{trn}$($v$)")
ax.plot(data[0], data[1], label=r"$s_{trn}^{g}(v)$")
ax.plot(data[0], data[2], label=r"$s_{trn}^{s}(v)$")
ax.set_xlim([0,500])
ax.set_xlabel("$v$ (cm$^{-1}$)", size="x-large")
ax.set_ylabel("s($v$) (cm)", size="x-large")
ax.legend()

ax = fig.add_subplot(3,2,5)
ax.plot(data[0], data[6], label="s$_{rot}$($v$)")
ax.plot(data[0], data[3], label="$s_{rot}^g(v)$")
ax.plot(data[0], data[4], label="$s_{rot}^s(v)$")
ax.set_xlim([0,1000])
ax.set_xlabel("$v$ (cm$^{-1}$)", size="x-large")
ax.set_ylabel("s($v$) (cm)", size="x-large")
ax.legend()

plt.legend()
plt.show()
