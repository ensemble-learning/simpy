""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import os
from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo

os.chdir("/home/tao/Documents/wag/alpha/opt")
a = Poscar("alpha1.vasp")
b = a.parser()
toXyz(b)
toGeo(b)

