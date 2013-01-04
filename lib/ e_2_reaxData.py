""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
from mytype import System, Molecule, Atom
from geo import Geo
from output_conf import toData
from output_conf import toGeo
from output_conf import toXyz

testfile = "../../debug/geo"
a = Geo(testfile)
b = a.parser()
b.assignAtomTypes()
toData(b)
toGeo(b)
toXyz(b)
