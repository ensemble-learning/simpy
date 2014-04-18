"""
Check the element parts of geo.
"""
import shutil
f = open("geo", "r")

lines = []

n = 0

for i in f:
    line = i
    if "HETATM" in i:
        ele = i[64:66]
        if len(ele.strip()) == 0:
            n += 1
            ele = i[13:15]
            line = i[:64] + ele + i[66:]
    lines.append(line)
f.close()

if n > 0:
    shutil.copy("geo", "geo.bak")
    print("Warning: The elements parts are missing in the GEO file"

o = open("geo", "w")
for i in lines:
    o.write(i)
o.close()

