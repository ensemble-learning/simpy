import os
import numpy as np
import matplotlib.pyplot as plt

x = []
y = []
yf_1 = []
yf_2 = []

f = open("rdf.xvg", "r")

for i in f:
    tokens = i.strip().split()
    if i.strip().startswith("#"):
        pass
    elif i.strip().startswith("@"):
        if "subtitle" in i:
            print tokens
            a1 = tokens[2]
            a2 = tokens[4]
    else:
        if len(tokens) == 2:
            x.append(float(tokens[0])*10)
            y.append(float(tokens[1]))
f.close()

#plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(x, y, lw=3)
#ax.plot(x, yf_1,ls="--")
ax.set_ylabel(r"$g_{Li-O}$ (r)", size=24)
ax.xaxis.grid(which="minor")
ax.xaxis.grid()

if os.path.exists("rdf_cn.xvg"):
    x = []
    y = []
    f = open("rdf_cn.xvg", "r")

    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("#"):
            pass
        elif i.strip().startswith("@"):
            pass
        else:
            if len(tokens) == 2:
                x.append(float(tokens[0])*10)
                y.append(float(tokens[1]))
    f.close()

    ax2 = ax.twinx()
    ax2.plot(x, y, lw=3, ls="--", color="black")
    ax2.set_ylabel(r"$\int g_{Li-O}$ (N)", size=24)
    ax2.set_ylim([0,15])

ax.set_xlabel("r ($\AA$)", size=24)
ax.set_xlim([0.5, 10.5])
plt.minorticks_on()
plt.tight_layout()
plt.savefig("figure.png", dpi=500)
plt.show()
