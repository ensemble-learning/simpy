"""
ITEM: ATOMS id type x y z vx vy vz
1 1 -49.9387 -60.8575 31.2613 0.00733834 -0.000396848 5.58796e-05
"""

f = open("dump.lammpstrj", "r")

for i in f: 
    if i.strip().startswith("ITEM: NUMBER OF ATOMS"):
        break

for i in f: 
    if i.strip().startswith("ITEM: BOX BOUNDS pp pp pp"):
        break
    else:
        natom = int(i.strip())

boundary = []
for i in f:
    if i.strip().startswith("ITEM: ATOMS"):
        break
    else:
        boundary.append([float(j) for j in i.strip().split()])

o = open("tmp.data", "w")
for i in f:
    tokens = i.strip().split()
    if len(tokens) > 5:
        n = int(tokens[0])
        id = int(tokens[1])
        x = float(tokens[2])
        y = float(tokens[3])
        z = float(tokens[4])
        charge = 0.0
        line = "%-10d%6d%8.4f%12.4f%12.4f%12.4f\n"%(n, id, charge, x, y, z)
        o.write(line)
o.close()
f.close()
        
