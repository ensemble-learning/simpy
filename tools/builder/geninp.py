""" read the pdb file and output to data (LAMMPS).
"""
import sys
import socket
import math
import os

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo
from cons import ATOMIC_MASS_CONSTANT as amc

def usage():
    print """python geninp pdbfile nmols density
    pdbfile: pdb file of single molecule
    nmols: number of molecules
    density: desired density
    """

def getBoxl(pdbfile, n, density):

    density = density * 1000
    a = Pdb(pdbfile)
    b = a.parser()
    b.getMass()
    mass = b.mass
    print "Molecule mass: ", mass
    total = mass * n
    vol = total*amc/density * 1e30 # here length is in A
    print "Volume: ", vol
    l = math.pow(vol, 1/3.0)
    print "box a, b, c =", l

    return l

def writeInp(pdbfile, n, l):
    from inp_tpl import INP
    lines = INP
    o = open("inp", "w")
    lines = lines.replace("%pdbfile%", pdbfile)
    lines = lines.replace("%n%", "%d"%n)
    lines = lines.replace("%l%", "%.4f"%l)
    o.write(lines)

def main():
    if len(sys.argv) < 4:
        usage()
    else:
        pdbfile = sys.argv[1]
        n = int(sys.argv[2])
        density = float(sys.argv[3])
        l = getBoxl(pdbfile, n, density)
        writeInp(pdbfile, n, l)
        os.system("packmol < inp")
        if os.path.exists("sim.pdb"):
            os.system("editconf -f sim.pdb -o sim.pdb -box %.4f"%(l/10))
            a = Pdb("sim.pdb")
            b = a.parser()
            b.assignAtomTypes()
            toReaxLammps(b)

if __name__ == "__main__":
    main()
