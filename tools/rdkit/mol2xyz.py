import sys

f = open('%s.mol'%sys.argv[1], 'r')
lines = f.readlines()
f.close()

tokens = lines[3].strip().split()
n_atoms = int(tokens[0])

coords = []
for i in range(len(lines)):
    tokens = lines[i].strip().split()
    if len(tokens) == 16:
        line = '%s %s %s %s\n'%(tokens[3], tokens[0], tokens[1], tokens[2])
        coords.append(line)

o = open('%s.xyz'%sys.argv[1], 'w')
o.write('%d\n\n'%len(coords))
for i in coords:
    o.write(i)
o.close()



