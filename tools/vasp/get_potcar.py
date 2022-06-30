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
elif socket.gethostname() == "fermion.wag.caltech.edu":
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
elif socket.gethostname() == "mu01":
    LIB = "/opt/sourcecoude/simpy/lib"
elif "stampede2" in socket.gethostname():
    LIB = "/home1/04076/tg833760/soft/simpy/lib"
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif "Tao-MBP-9" in socket.gethostname():
    LIB = "/Users/tao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb

"""
POT = {"N":"N", "O":"O", "H":"H", "C":"C", "Li":"Li", "S":"S", "Ti":"Ti", "P":"P", 
       "Ca":"Ca_pv", "Al":"Al", "Cu":"Cu", "Na":"Na", "Cl":"Cl", "Cr":"Cr", "Ga":"Ga",
        "Br":"Br", "D": "H", "Si": "Si", "Ni": "Ni", "Pt":"Pt_pv", "Co":"Co", "Cr":"Cr", 
       "I":"I", "K":"K", "F":"F", "W":"W", "Au":"Au", "Cs":"Cs_sv", "Mg":"Mg",
       "Ag":"Ag", "Se":"Se", "B":"B", "He":"He", "Ar":"Ar", "Xe":"Xe", "Kr":"Kr", 
       "Mo": "Mo", "Fe": "Fe", "As": "As", "Ge": "Ge", "Sc": "Sc", "Zr": "Zr",
       "Y": "Y_sv", "Zn": "Zn", "Cd": "Cd", "V": "V", "Mn": "Mn", "Co": "Co",
       "Nb": "Nb_sv", "Tc": "Tc", "Ru": "Ru", "Rh": "Rh_pv", "Pd": "Pd", "Ne": "Ne", 
       "Bi": "Bi", "Ba": "Ba_sv", "Hf": "Hf", "In":"In", "Ir":"Ir", "Lu":"Lu", "Os":"Os",
       "Pb":"Pb", "Te":"Te", "Re":"Re", "Sb":"Sb", "Sn":"Sn", "Sr":"Sr_sv", "Ta": "Ta", "Te":"Te",
       "Tl":"Tl", "Po":"Po", "Be":"Be"}
"""
# for MPL
POT = {"N":"N", "O":"O", "H":"H", "C":"C", "Li":"Li_sv", "S":"S", "Ti":"Ti_pv", "P":"P", 
       "Ca":"Ca_sv", "Al":"Al", "Cu":"Cu_pv", "Na":"Na_pv", "Cl":"Cl", "Ga":"Ga_d",
        "Br":"Br", "D": "H", "Si": "Si", "Ni": "Ni_pv", "Pt":"Pt_pv", "Co":"Co", "Cr":"Cr_pv", 
       "I":"I", "K":"K_pv", "F":"F", "W":"W_pv", "Au":"Au", "Cs":"Cs_sv", "Mg":"Mg_pv",
       "Ag":"Ag", "Se":"Se", "B":"B", "He":"He", "Ar":"Ar", "Xe":"Xe", "Kr":"Kr", 
       "Mo": "Mo_pv", "Fe": "Fe_pv", "As": "As", "Ge": "Ge_d", "Sc": "Sc_sv", "Zr": "Zr_sv",
       "Y": "Y_sv", "Zn": "Zn", "Cd": "Cd", "V": "V_pv", "Mn": "Mn_pv", "Co": "Co",
       "Nb": "Nb_pv", "Tc": "Tc_pv", "Ru": "Ru_pv", "Rh": "Rh", "Pd": "Pd", "Ne": "Ne", 
       "Bi": "Bi", "Ba": "Ba_sv", "Hf": "Hf_pv", "In":"In_d", "Ir":"Ir", "Lu":"Lu_3", "Os":"Os_pv",
       "Pb":"Pb_d", "Te":"Te", "Re":"Re_pv", "Sb":"Sb", "Sn":"Sn_d", "Sr":"Sr_sv", "Ta": "Ta_pv", "Te":"Te",
       "Tl":"Tl_d", "Po":"Po", "Be":"Be_sv", "Ac": "Ac", "Ce":"Ce", "Dy":"Dy_3", "Er": "Er_3",
       "Gd":"Gd", "Hg":"Hg", "Ho":"Ho_3", "La":"La", "Nd": "Nd_3", "Np": "Np", "Pa": "Pa", 
       "Pm": "Pm_3", "Pr": "Pr_3", "Pu": "Pu", "Rb": "Rb_sv", "Sm": "Sm_3", "Tb": "Tb_3", 
       "Th": "Th", "Tm": "Tm_3", "U": "U", "Yb": "Yb_2", "Eu":"Eu"}
 
if socket.gethostname() == "cluster.hpc.org":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "tao-laptop":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "atom.wag.caltech.edu":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "ion.wag.caltech.edu":
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "fermion.wag.caltech.edu":
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
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    POT_DATA_BASE = "/home/tao/data/soft/vasp/vasp.5.3.5/potcar/potpaw_PBE"
elif socket.gethostname() == "mu05":
    POT_DATA_BASE = "/home/chengtao/soft/vasp/potpaw_PBE"
elif socket.gethostname() == "mu01":
    POT_DATA_BASE = "/opt/software/vasp.5.4.4/potpaw_PBE"
elif "stampede2" in socket.gethostname():
    POT_DATA_BASE = "/home1/04076/tg833760/soft/vasp/potcar/potpaw_PBE"
elif "node" in socket.gethostname():
    POT_DATA_BASE = "/project/source/VASP/vasp.5.3.5/potcar/potpaw_PBE"
elif "Tao-MBP-9" in socket.gethostname():
    POT_DATA_BASE = "/Users/tao/soft/vasp/vasp.5.4.4/potpaw_PBE"

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
