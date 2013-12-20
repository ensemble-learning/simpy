f = open("geo_end.xyz", "r")
o = open("lammps.xyz", "w")

for i in f:
    tokens = i.strip().split()
    #@note: not a strick way to get the coords
    if len(tokens) == 8:
        line = " ".join(tokens[:4]) + "\n"
    else:
        line = i
    o.write(line)
o.close()
f.close()
