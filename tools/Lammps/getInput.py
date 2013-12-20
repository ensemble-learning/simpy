""" This code is used to generate reaxFF lammps_input file
according to lammps.data (also generated from simpy
"""
import sys
import os
import socket
import argparse

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simupy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from template import *
from ffield import Ffield

#C H O N
#FF = {"C": 1, "H": 2, "O": 3, "N": 4, "Ca":4, "Al":6}
FF = {}
MASS = {12.011:"C", 14.007: "N", 15.994:"O", 1.0079:"H", 40.078:"Ca",\
        26.982:"Al", 28.086:"Si", 35.453:"Cl", 47.867:"Ti",  6.941:"Li",
        30.974:"P", 32.065:"S", 72.64:"Ge", 10.811:"B"}

def usage():
    print """python genInput.py type
    type: NVT, MIN
    """

def getElements(args):
    assert os.path.exists("ffield")
    if args.lg:
        a = Ffield("ffield", 1 )
    else:
        a = Ffield("ffield", 0 )
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

def main(args):
    """ generate the lammps input file
    """
    getElements(args)
    m = parseData()
    ty = []
    elem = []
    for i in m:
        #i = int(i*1000)/1000.0
        at = MASS[i]
        elem.append(at)
        ff = FF[at]
        ty.append("%d"%ff)

    rtype = args.type
    lg = 0
    if args.lg:
        lg = 1

    if rtype == "NVT":
        lines = NVT
    elif rtype == "MIN":
        lines = MIN
    elif rtype == "NPT":
        lines = NPT
    elif rtype == "MIN_CELL":
        lines = MIN_CELL

    print "processing %s simulation......"%rtype
    
    if lg:
        lines = lines.replace("%reax_potential%", "reax/c NULL lgvdw yes")
    else:
        lines = lines.replace("%reax_potential%", "reax/c NULL")
    
    if args.lammps2012:
        lines = lines.replace("%ffield_atoms%", " ".join(ty))
    else:
        lines = lines.replace("%ffield_atoms%", " ".join(elem))

    lines = lines.replace("%elements%", " ".join(elem))

    o = open("lammps_input", "w")
    o.write(lines)
    o.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("type", default="MIN", nargs="?", help="simulation type")
    parser.add_argument("-lg", action="store_true", help="using lg type ffield")
    parser.add_argument("-lammps2012", action="store_true", help="using lg type ffield")
    #parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    #parser.add_argument("-b", nargs=2, type=int, help="get the bond distance between a1, a2, a3")
    #parser.add_argument("-a", nargs=3, type=int,help="get the angle of a1-a2-a3")
    #parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")
    args = parser.parse_args()
    main(args)


