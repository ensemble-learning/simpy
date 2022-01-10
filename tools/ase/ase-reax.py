import os, sys
from ase.io import read, write
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS, LBFGS, GPMin, FIRE
from numpy import sort, unique
from mendeleev import element

fname = sys.argv[1]

atoms = read(fname)
if not all(atoms.get_pbc()):
    atoms.set_pbc([True, True, True])
    atoms.set_cell([50,50,50])

elems = unique(sort(atoms.get_atomic_numbers()))

symbols = []
for i in elems:
    symbols.append((element(int(i)).symbol))

# need to pay special attention to the order of the atoms
parameters = {
        'units': 'real',
        'atom_style': 'charge',
        'pair_style': 'reax/c control.reaxc',
        'pair_coeff': ['* * ffield %s'%' '.join(symbols)],
        'min_style': 'cg',
        'minimize': '0 1.0e-8 1000 1000',
        'fix':['QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c'],
        }

files = ['control.reaxc', 'ffield']
lammps = LAMMPS(parameters=parameters,
                        files=files, tmp_dir='lammps-run',
                        )

atoms.calc = lammps
#opt = LBFGS(atoms)
#opt.run(fmax=0.02, steps=1000)

pe = atoms.get_potential_energy()
write('reax-out.pdb', atoms)


