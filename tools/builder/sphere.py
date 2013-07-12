""" Cut a sphere from a pdb input file
@testcase: example 01
@log:
"""

import sys, os
import socket
import math
import ConfigParser, string
import copy

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo, toPdb
from cons import ATOMIC_MASS_CONSTANT as amc
from utilities import get_dist

# max atoms allowed
MAX_ATOMS = 99999999

class SphereCut():
    """Cut a sphere from a pdb file
    """
    def __init__(self,):
        self.center = [0.0, 0.0, 0.0]
        self.radius = 0.0
        self.pdbfile = ''
        self.atoms = []
        self.natoms = []
        self.output = "tmp.pdb"

def read_inp(control):
    """Read the input file
    """
    print "Reading the control file ...."
    cf = ConfigParser.ConfigParser()
    cf.read("inp")

    s = cf.sections()
    if "SPHERE" in s:
        print "    Initialize sphere cut ..."
        o = cf.options("SPHERE")
        if "pdbfile" in o:
            control.pdbfile = cf.get("SPHERE", "pdbfile").strip()
            print "        model file is %s"%control.pdbfile
        if "center" in o:
            tokens = cf.get("SPHERE", "center").strip().split()
            control.center[0] = float(tokens[0])
            control.center[1] = float(tokens[1])
            control.center[2] = float(tokens[2])
            print "        center of sphere is %.3f %.3f %.3f"%(control.center[0],\
            control.center[1], control.center[2])
        if "radius" in o:
            control.radius = float(cf.get("SPHERE", "radius").strip())
            print "        Radius is %.3f"%control.radius
        if "atoms" in o:
            control.atoms = cf.get("SPHERE", "atoms").strip().split()
        if "natoms" in o:
            tokens = cf.get("SPHERE", "natoms").strip().split()
            control.natoms = [int(i) for i in tokens]
        if "natoms" in o:
            tokens = cf.get("SPHERE", "natoms").strip().split()
            control.natoms = [int(i) for i in tokens]
        if "outfile" in o:
            control.output = cf.get("SPHERE", "outfile").strip()

def read_pdb(control):
    print "    Reading pdb file ..."
    a = Pdb(control.pdbfile)
    b = a.parser()
    b.assignAtomTypes2()
    print "        cell parameters:"
    print "            a = %8.3f: "%b.pbc[0]
    print "            a = %8.3f: "%b.pbc[1]
    print "            a = %8.3f: "%b.pbc[2]
    return b

def sphere_cut(control, b):
    """Cut the sphere and balance the charge.
    Here, we use the ratio from input as criteria. We eliminate the extra 
    atoms according to the distances to the center. (So we first sort the 
    atoms according to distance).
    """
    
    # coarse cut (only according to radius
    print "    Cutting the sphere ..."
    natoms = [0]*len(control.atoms)
    ndx1 = []
    ndx2 = []
    for i in range(len(control.atoms)):
        ndx1.append([])
        ndx2.append([])

    counter = 0
    for i in b.atoms:
        dist = get_dist(i.x, control.center)
        if dist < control.radius:
            n = control.atoms.index(i.name)
            natoms[n] += 1
            ndx1[n].append("%09d_%08d"%(dist*1000, counter))
        counter += 1

    # print information of coarse cut 
    print "        After Coarse cut we get:"
    for i in range(len(control.atoms)):
        print "            %-6s = %8d"%(control.atoms[i], len(ndx1[i]))

    # sort the atoms
    for i in ndx1:
        i.sort()

    nmin = MAX_ATOMS
    for i in range(len(control.natoms)):
        n = int(natoms[i]/control.natoms[i])
        if n < nmin:
            nmin = n

    for i in range(len(control.atoms)):
        n = nmin * control.natoms[i]
        for j in range(n):
            ndx2[i].append(int(ndx1[i][j].split("_")[-1]))

    print "        After balancing charge we get:"
    for i in range(len(control.atoms)):
        print "            %-6s = %8d"%(control.atoms[i], len(ndx2[i]))
    return ndx2
    
def output_pdb(control, b, ndx):
    """Output the atoms in ndx to pdb file
    """
    c = copy.deepcopy(b)
    c.atoms = []
    for i in range(len(ndx)):
        for j in range(len(ndx[i])):
            c.atoms.append(b.atoms[ndx[i][j]])
    toPdb(c, control.output)

def main():

    control = SphereCut()
    # read input
    read_inp(control)
    # read pdb
    b = read_pdb(control)
    # sphere cut
    ndx = sphere_cut(control, b)
    # output
    output_pdb(control, b, ndx)

if __name__ == "__main__":
    main()

