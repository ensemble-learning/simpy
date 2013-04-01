"""modify the description section
"""

import shutil

lines = []
f = open("geo", "r")

for i in f:
    if i.strip().startswith("DESCR"):
        tokens = i.strip().split()
        tokens[1] = "r_hdn_" + tokens[1]
        line = " ".join(tokens) + "\n"
    else:
        line = i
    lines.append(line)
f.close()
shutil.copy("geo", "geo.bak")

o = open("geo", "w")
for i in lines:
    o.write(i)
o.close()

