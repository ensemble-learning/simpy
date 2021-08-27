from ase.io import read, write
A2B = 1.889726126
atoms = read('POSCAR')

# get info
ele = atoms.get_chemical_symbols()
pos = atoms.get_positions()
cell = atoms.get_cell()
cell = cell.transpose() # beware that the sequence is different from vasp

# write lattice
o = open('lattice', 'w')
o.write('lattice \ \n')
for i in cell:
    for j in i:
        o.write('%18.9f '%(j*A2B))
    o.write('\ \n')
o.close()

# write coords
o = open('ionpos', 'w')
for n, i in enumerate(ele):
    o.write('ion    %2s '%i)
    for j in pos[n]:
        o.write('%12.6f '%(j*A2B))
    o.write(' 0\n') # default fix all the atoms, assuming single point calculation
o.close()

