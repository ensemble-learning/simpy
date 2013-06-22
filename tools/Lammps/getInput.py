""" This code is used to generate reaxFF lammps_input file
according to lammps.data (also generated from simpy
"""
import sys
import os
import socket

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simupy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"

sys.path.insert(0 , LIB)

from template import *
from ffield import Ffield

#C H O N
#FF = {"C": 1, "H": 2, "O": 3, "N": 4, "Ca":4, "Al":6}
FF = {}
MASS = {12.011:"C", 14.007: "N", 15.994:"O", 1.0079:"H", 40.078:"Ca",\
        26.982:"Al"}

def usage():
    print """python genInput.py type
    type: NVT, MIN
    """

def getElements():
    assert os.path.exists("ffield")
    a = Ffield("ffield")
    counter = 1
    for i in a.elements:
        FF[i] = counter
        counter += 1

def parseData(fname="lammps.data"):
    """parse the data file to get the mass infor
    """
    m = []
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("Masses"):
            break
    for i in f:
        if len(i.strip()) == 0:
            pass
        elif i.strip().startswith("Atoms"):
            break
        else:
            tokens = i.strip().split()
            amass = float(tokens[1])
            m.append(amass)
    f.close()
    return m

def main(rtype="MIN"):
    """ generate the lammps input file
    """
    getElements()
    m = parseData()
    ty = []
    for i in m:
        at = MASS[i]
        ff = FF[at]
        ty.append("%d"%ff)

    if rtype == "NVT":
        lines = NVT
    elif rtype == "MIN":
        lines = MIN
    elif rtype == "NPT":
        lines = NPT
    elif rtype == "MIN_CELL":
        lines = MIN_CELL

    print "processing %s simulation......"%rtype
    lines = lines.replace("%ffield_atoms%", " ".join(ty))
    o = open("lammps_input", "w")
    o.write(lines)
    o.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
        print "Waring: using default input type (MIN)"
        main()
    else:
        main(sys.argv[1])
