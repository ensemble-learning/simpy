import time
import os

rdfs = []
f = open('rdf.dat', 'r')
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 2:
        rdfs.append([tokens[0], tokens[1]])
f.close()

for i in rdfs:
    o = open('inp-ndx.dat', 'w')
    o.write('a %s\n'%i[0])
    o.write('a %s\n'%i[1])
    o.write('q\n')
    o.write('\n')
    o.close()

    time.sleep(0.1)
    os.system('cat inp-ndx.dat | gmx make_ndx -f input.pdb')

    o = open('inp-rdf.dat', 'w')
    o.write('%s\n'%i[0])
    o.write('%s\n'%i[1])
    o.write('\n')
    o.write('\n')
    o.close()

    time.sleep(0.1)
    os.system('cat inp-rdf.dat | gmx rdf -f output.pdb -n index.ndx -o %s-%s -cn %s-%s-cn'%(i[0], i[1], i[0], i[1]))

