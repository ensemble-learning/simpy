from ase.build import molecule
from ase.collections import g2
from ase.io import write

for i in g2.names:
    atoms = molecule(i)
    write(i+'.pdb', atoms)
    
