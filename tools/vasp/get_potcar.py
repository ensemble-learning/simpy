import sys
import os
import socket

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb

POT = {"N":"N", "O":"O", "H":"H", "C":"C", "Li":"Li", "S":"S", "Ti":"Ti", "P":"P", 
       "Ca":"Ca_pv", "Al":"Al", "Cu":"Cu", "Na":"Na", "Cl":"Cl", "Cr":"Cr", "Ga":"Ga",
        "Br":"Br", "D": "H", "Si": "Si", "Ni": "Ni", "Pt":"Pt"}
POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"

o = open("POTCAR", "w")
a = Poscar("POSCAR")
for i in a.atomtypes:
    print i
    pot = os.path.join(POT_DATA_BASE, POT[i])
    pot = os.path.join(pot, "POTCAR")
    f = open(pot, "r")
    for j in f:
        o.write(j)
    f.close()
o.close()
    
    
    


