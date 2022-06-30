from ase.io import read, write
from ase.build import sort

atoms = read('pt.pdb')
com = atoms.get_center_of_mass()
cell = atoms.get_cell()
dx = []
for i in range(3):
    dx.append(-com[i] + cell[i][i]/2)

atoms.translate(dx)
atoms = sort(atoms)

write('input.pdb', atoms)
write('POSCAR', atoms, vasp5=True, direct=True)
