from ase.io import read, write
import sys
from ase.build import sort

fname = sys.argv[1]
catom = int(sys.argv[2])

atoms = read(fname)
sp = atoms.get_scaled_positions()
dx = [0.5, 0.5, 0.5] - sp[catom-1]
for n, i in enumerate(sp):
    sp[n] = sp[n] + dx
atoms.set_scaled_positions(sp)
sp = atoms.get_scaled_positions()
atoms.wrap()
atoms = sort(atoms)
write('out.vasp', atoms, vasp5=True, direct=True)
