import sys

f = open('%s.mol'%sys.argv[1], 'r')
lines = f.readlines()
f.close()

tokens = lines[3].strip().split()
n_atoms = int(tokens[0])

coords = []
for i in range(4, 4+n_atoms):
    tokens = lines[i].strip().split()
    line = '%s %s %s %s\n'%(tokens[3], tokens[0], tokens[1], tokens[2])
    coords.append(line)

o = open('%s.xyz'%sys.argv[1], 'w')
o.write('%d\n\n'%n_atoms)
for i in coords:
    o.write(i)
o.close()



