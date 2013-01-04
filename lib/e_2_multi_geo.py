""" parse the geo file with multi configuration into
seperated files
"""
import os
from utilities import parseBlock
from mytype import System, Molecule, Atom
from geo import Geo
from output_conf import toGeo
from output_conf import toXyz

os.chdir("/home/tao/Documents/debug/geofile")
parseBlock("geo", 1)
for i in range(204):
    fname = "out%03d"%i 
    a = Geo(fname)
    b = a.parser()
    b.assignAtomTypes()
    toGeo(b, b.name+'.geo')
    toXyz(b, b.name+'.xyz')

