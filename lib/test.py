from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toPoscar

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

toPoscar(b)


