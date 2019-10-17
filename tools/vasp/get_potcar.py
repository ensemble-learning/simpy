#!/usr/bin/env python

import sys
import os
import socket

LIB = ''

print(socket.gethostname())

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
elif "comet" in socket.gethostname():
    LIB = "/home/tcheng/soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "zwicky":
    LIB = "/home/tcheng/Soft/simpy/lib"
elif "onyx" in socket.gethostname():
    LIB = "/p/home/taocheng/src/simpy/lib"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    LIB = "/home/tao/soft/simpy/lib"
elif socket.gethostname() == "mu05":
    LIB = "/home/chengtao/soft/simpy/lib"
elif "stampede2" in socket.gethostname():
    LIB = "/home1/04076/tg833760/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb

POT = {"N":"N", "O":"O", "H":"H", "C":"C", "Li":"Li", "S":"S", "Ti":"Ti", "P":"P", 
       "Ca":"Ca_pv", "Al":"Al", "Cu":"Cu", "Na":"Na", "Cl":"Cl", "Cr":"Cr", "Ga":"Ga",
        "Br":"Br", "D": "H", "Si": "Si", "Ni": "Ni", "Pt":"Pt_pv", "Co":"Co", "Cr":"Cr", 
       "I":"I", "K":"K_pv", "F":"F", "W":"W", "Au":"Au", "Cs":"Cs_sv", "Mg":"Mg",
       "Ag":"Ag", "Se":"Se", "B":"B", "He":"He", "Ar":"Ar", "Xe":"Xe", "Kr":"Kr", 
       "Mo": "Mo", "Fe": "Fe", "As": "As", "Ge": "Ge", "Sc": "Sc", "Zr": "Zr",
       "Y": "Y_sv", "Zn": "Zn", "Cd": "Cd", "V": "V", "Mn": "Mn", "Co": "Co",
       "Nb": "Nb_sv", "Tc": "Tc", "Ru": "Ru", "Rh": "Rh", "Pd": "Pd", "Ne": "Ne", 
       "Bi": "Bi", "Ba": "Ba_sv", "Hf": "Hf", "In":"In", "Ir":"Ir", "Lu":"Lu", "Os":"Os",
       "Pb":"Pb", "Te":"Te", "Re":"Re", "Sb":"Sb", "Sn":"Sn", "Sr":"Sr_sv", "Ta": "Ta", "Te":"Te",
       "Tl":"Tl", "Po":"Po"}

if socket.gethostname() == "cluster.hpc.org":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "tao-laptop":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "atom.wag.caltech.edu":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "ion.wag.caltech.edu":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "giant12":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif "comet" in socket.gethostname():
    POT_DATA_BASE = "/home/tcheng/soft/vasp.5.3/potpaw_PBE"
elif socket.gethostname() == "tao-Precision-Tower-3420":
    POT_DATA_BASE = "/home/tao/Soft/vasp/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "zwicky":
    POT_DATA_BASE = "/home/tcheng/Soft/potpaw_PBE"
elif "onyx" in socket.gethostname():
    POT_DATA_BASE = "/p/home/taocheng/src/vasp/vasp.5.3/potcar/potpaw_PBE"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    POT_DATA_BASE = "/home/tao/src/vasp/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    POT_DATA_BASE = "/home/tao/src/vasp/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    POT_DATA_BASE = "/home/tao/data/soft/vasp/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "mu05":
    POT_DATA_BASE = "/home/chengtao/soft/vasp/potpaw_PBE"
elif "stampede2" in socket.gethostname():
    POT_DATA_BASE = "/home1/04076/tg833760/soft/vasp/potcar/potpaw_PBE"

o = open("POTCAR", "w")
a = Poscar("POSCAR")
for i in a.atomtypes:
    print(i)
    pot = os.path.join(POT_DATA_BASE, POT[i])
    pot = os.path.join(pot, "POTCAR")
    f = open(pot, "r")
    for j in f:
        if j.strip().startswith("VR"):
            print(j)
        o.write(j)
    f.close()
o.close()
    
    
    


