from lib.mytype import System, Molecule, Atom
from lib.geo import Geo
from lib.output import toData
from lib.output import toGeo
from lib.output import toXyz

testfile = "../debug/geo"
a = Geo(testfile)
b = a.parser()
b.assignAtomTypes()
toData(b)
toGeo(b)
toXyz(b)

