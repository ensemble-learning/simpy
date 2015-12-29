import sys, socket
import copy
import os
import math

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
elif socket.gethostname() == "giant3":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant1":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from data import ReaxData
from output_conf import toPdb, toReaxLammps

CN = """[GLOBAL]
r_ref_0 = 2.788
r_ref_1 = 3.944
r_ref_2 = 4.830

n_atoms = %natoms%
# cut_off_func = 1(cut-off) 2(cn function)
cut_off_func = 1
r_cut_range = 1.8
n_cut_off = 4
"""

SAS = """UFF.atoms
out.music
1.400
5000     
%pbc%
0.59
0.0
0.04
cn.dat
nlist.dat



! file containing the atom types + diameters
! file containing the coordinates
! probe size in A
! number of insertions
! length of unitcell
! crystal density in g / cm3
! scale factor cut off
! surface atom cut off
! coordinatin files
! neighbor list files
"""

datafile = "lammps.data"
a = ReaxData(datafile)
b = a.parser()
b.assignAtomTypes()

natoms = len(b.atoms)
pbc = b.pbc

o = open("cn.inp", "w")
lines = CN
lines = lines.replace("%natoms%", "%d"%natoms)
o.write(lines)
o.close()

o = open("sas.inp", "w")
lines = SAS
tokens = " ".join(["%.4f"%i for i in pbc[:3]])
lines = lines.replace("%pbc%", tokens)
o.write(lines)
o.close()

o = open("run.sh", "w")

o.write("""
lammps_debug -in lammps_input_debug
sleep 1
python ~/soft/simpy/tools/sas/cn.py
sleep 1
sas < sas.inp
""")
o.close()

os.system("chmod +x run.sh")
