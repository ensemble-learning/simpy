from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb, toReaxLammps

import sys

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = "CONTCAR"
a = Poscar(fname)
b = a.parser()
b.assignAtomTypes()
b.assignEleTypes()
print b.getVol()

toXyz(b)
toGeo(b, "geo")
toPdb(b, "sim.pdb")
toReaxLammps(b, "lammps.data")

