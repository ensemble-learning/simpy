""" parse the geo file with multi configuration into
seperated files
"""
import os
from utilities import parseBlock
from mytype import System, Molecule, Atom
from geo import Geo
from output_conf import toGeo, toXyz, toReaxLammps, toPdb

fname = "geo"
a = Geo(fname)
b = a.parser()
b.assignAtomTypes()
toGeo(b, b.name+'.geo')
toXyz(b, b.name+'.xyz')
toReaxLammps(b)
toPdb(b)

