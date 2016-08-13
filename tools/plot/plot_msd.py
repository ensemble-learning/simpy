import numpy as np
import matplotlib.pyplot as plt

x = []
y = []
yf_1 = []
yf_2 = []

f = open("msd.xvg", "r")

for i in f:
    tokens = i.strip().split()
    if i.strip().startswith("#"):
        if "cm^2/s" in i:
            D = float(tokens[4])
    elif i.strip().startswith("@"):
        pass
    else:
        if len(tokens) == 2:
            x.append(float(tokens[0]))
            y.append(float(tokens[1]))
f.close()

for i in range(len(x)):
    v1 = 6*D*x[i]/1000
    yf_1.append(v1)

plt.figure(figsize=(10,5))
ax = plt.subplot(121)
ax.plot(x, y, lw=3)
#ax.plot(x, yf_1,ls="--")
ax.set_xlabel("T (ps)", size=24)
ax.set_ylabel("MSD (nm$^2$)", size=24)
ax.set_title("MSD")

ax = plt.subplot(122)
ax.plot(x, y, lw=3)
ax.plot(x, yf_1,ls="--", label="~t", lw=3)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("T (ps)", size=24)
#ax.set_ylabel("MSD (nm$^2$)", size=24)
ax.set_title("MSD (log vs log)")
ax.legend(loc=2, fontsize=18)
plt.tight_layout()
plt.savefig("figure.png", dpi=500)
plt.show()
