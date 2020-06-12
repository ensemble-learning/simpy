from math import pi
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
#from deepmd.calculator import DP
from ase.build import molecule
from ase.constraints import FixInternals

atoms = read('input.vasp')
for i in atoms:
    for n in range(3):
        i.position[n] = i.position[n] + 5
atoms.wrap()

write('POSCAR', atoms, vasp5=True)

