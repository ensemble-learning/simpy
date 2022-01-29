#!/usr/bin/env python

import os, sys
from ase.io import read, write
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS, LBFGS, GPMin, FIRE
from numpy import sort, unique
from mendeleev import element
import shutil

ff = '/home/tcheng/share/reaxff/ffield.reax.cho'
control = '/home/tcheng/share/reaxff/control.reaxc'

fname = sys.argv[1]

shutil.copy(ff, 'ffield')
shutil.copy(control, 'control.reaxc')

atoms = read(fname)
com = atoms.get_center_of_mass()

if not all(atoms.get_pbc()):
    dx = []
    atoms.set_pbc([True, True, True])
    atoms.set_cell([50,50,50])
    for i in range(3):
        dx.append(-com[i] + 25)
    atoms.translate(dx)

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
o = open('reax-ener.dat', 'w')
o.write('%.6f eV\n'%pe)
o.close()
write('reax-out.pdb', atoms)


