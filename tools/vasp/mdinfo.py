steps = []
pote = []
potf = []

f = open("OSZICAR", "r")

for i in f:
    if "=" in i:
        tokens = i.strip().split()
        steps.append(int(tokens[0]))
        pote.append(float(tokens[8]))

f.close()

o = open("potential", "w")
for i in pote:
    o.write("%.4f\n"%i)
o.close()

import matplotlib.pyplot as plt

plt.plot(pote)
plt.show()

