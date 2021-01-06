from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.neighborlist import natural_cutoffs

nt = [0]*32
atoms = read('run-01.pdb')
nl = NeighborList(natural_cutoffs(atoms, 1.35))
nl.update(atoms)

o = open('nl.dat', 'w')
for i in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    nt[len(indices)] += 1
    o.write('%d\n'%len(indices))

