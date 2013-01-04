""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
from mytype import System, Molecule, Atom
from data import ReaxData
from output_conf import toPdb

testfile = "e_2_data.data"
a = ReaxData(testfile)
b = a.parser()
b.assignAtomTypes()
toPdb(b)
