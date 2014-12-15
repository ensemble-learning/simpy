atoms = []

f = open("mol.atom", "r")

for i in f:
    if i.strip().startswith("#"):
        pass
    else:
        tokens = i.strip().split()
        if len(tokens) == 3:
            molname = tokens[1]
            id = int(tokens[0])
            if molname == "N2":
                atoms.append(id)
f.close()

o = open("n2.atoms", "w")
for i in atoms:
    o.write("%d\n"%i)
o.close()
