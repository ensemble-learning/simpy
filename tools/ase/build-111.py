import os
from ase.build import molecule
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.io import write

symbol = 'Pt'

folder = '01-slab'
if not os.path.exists(folder):
    os.mkdir(folder)
os.chdir(folder)
slab = fcc111(symbol, size=(3,3,4), orthogonal=False)
slab.center(vacuum=15.0, axis=2)
c = FixAtoms(indices=[atom.index for atom in slab if atom.symbol == symbol])
slab.set_constraint(c)
write('POSCAR', slab, vasp5=True, direct=True)
os.chdir('..')

folder = '02-h'
if not os.path.exists(folder):
    os.mkdir(folder)
os.chdir(folder)
slab = fcc111(symbol, size=(3,3,4), orthogonal=False)
add_adsorbate(slab, 'H', 1.0, 'fcc')
slab.center(vacuum=15.0, axis=2)
c = FixAtoms(indices=[atom.index for atom in slab if atom.symbol == symbol])
slab.set_constraint(c)
write('POSCAR', slab, vasp5=True, direct=True)
os.chdir('..')

folder = '03-h2o'
if not os.path.exists(folder):
    os.mkdir(folder)
os.chdir(folder)
slab = fcc111(symbol, size=(3,3,4), orthogonal=False)
water = molecule('H2O')
water.rotate(90, 'y')
add_adsorbate(slab, water, 2.0, 'fcc')
slab.center(vacuum=15.0, axis=2)
c = FixAtoms(indices=[atom.index for atom in slab if atom.symbol == symbol])
slab.set_constraint(c)
write('POSCAR', slab, vasp5=True, direct=True)
os.chdir('..')

"""
print(slab.get_tags())
"""
