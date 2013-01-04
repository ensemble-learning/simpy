from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb

a = Poscar("CONTCAR")
b = a.parser()
print b.getVol()

toXyz(b)
toGeo(b, "geo")
toPdb(b, "sim.pdb")

