from ase.build import bulk
from ase.io import write

for i in ['Pt', 'Ir', 'Pd', 'Rh']:
    atoms = bulk(i, 'fcc', cubic=True)
    write('%s.vasp'%i, atoms, vasp5=True)
    
