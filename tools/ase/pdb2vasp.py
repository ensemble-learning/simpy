import sys
from ase.io import read, write

if len(sys.argv) > 1:
    pdbfile = sys.argv[1]
    atoms = read(pdbfile)
    vaspfile = pdbfile.split('.')[0] + '.vasp'
    write(vaspfile, atoms, vasp5=True)

