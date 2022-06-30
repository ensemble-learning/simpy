atoms = []
new_atoms = []

f = open('atoms.0', 'r')
for i in f:
    tokens = i.strip()
    if len(tokens) > 0:
        atoms.append(tokens)
f.close()

h_atoms = []
n_no_h = 0
ndx =  {}
for n, i in enumerate(atoms):
    if i == 'H':
        h_atoms.append(str(n+1))
    else:
        ndx[str(n+1)] = str(n_no_h+1)
        new_atoms.append(i)
        n_no_h += 1

new_bonds = []
f = open('bonds.0', 'r')
for i in f:
    tokens = i.strip().split()
    if not ((tokens[0] in h_atoms) or (tokens[1] in h_atoms)):
        new_bonds.append([ndx[tokens[0]], ndx[tokens[1]], tokens[2]])
f.close()

o = open('atoms', 'w')
for i in new_atoms:
    o.write(i+'\n')
o.close()

o = open('bonds', 'w')
for i in new_bonds:
    o.write(' '.join(i)+'\n')
o.close()
    
