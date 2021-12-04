from ase.io import read, write
atoms = read('out.pdb')
com = atoms.get_center_of_mass()
cell = atoms.get_cell()
dx = []
for i in range(3):
    dx.append(-com[i] + cell[i][i]/2)

atoms.translate(dx)
write('pt-nw-new.pdb', atoms)
