import sys

infile = sys.argv[1] + ".gjf"
outfile = sys.argv[1] + ".xyz"

"""

%chk=D:\Dropbox\e-n-020.chk
# hf/3-21g geom=connectivity

CONTCAR

0 1
"""


f = open(infile, "r")

for i in f:
    if i.strip().startswith("%"):
        pass
    elif i.strip().startswith("#"):
        pass
    else:
        tokens = i.strip()
        if len(tokens) == 0:
            break

for i in f:
    tokens = i.strip()
    if len(tokens) == 0:
        break

atoms = []
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 4:
        atoms.append(tokens)
    if len(tokens) == 0:
        break
f.close()

o = open(outfile, "w")
o.write("%d\n"%(len(atoms)))
o.write("\n")

print atoms

for i in atoms:
    o.write("%6s%12.4f%12.4f%12.4f\n"%(i[0],
            float(i[1]), float(i[2]), float(i[3])))
o.close()
        

