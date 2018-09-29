n2_atoms = []
f = open("n2.atoms", "r")
for i in f:
    tokens = i.strip()
    n2_atoms.append(int(tokens))
f.close()

f = open("config.out", "r")
o = open("n2.lammpstrj", "w")

o .write("""ITEM: TIMESTEP
300000
ITEM: NUMBER OF ATOMS
%d
ITEM: BOX BOUNDS xy xz yz pp pp ff
-71.9148 160.692 -71.9148
0 126.748 0
0 1480.6 0
ITEM: ATOMS id type q x y z vx vy vz
"""%len(n2_atoms))

n = 0
for i in f:
    if i.startswith("#"):
        pass
    else:
        tokens = i.strip().split()
        if len(tokens) == 10:
            id = int(tokens[0])
            if n < len(n2_atoms):
                if id == n2_atoms[n]:
                    n += 1
                    line = tokens[0] + " 2 " + " ".join(tokens[3:]) + '\n'
                    o.write(line)
o.close()
f.close()
    
    

