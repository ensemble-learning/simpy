""" read the pdb file and output to data (LAMMPS).
"""
import sys
import socket

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo

def test():
    testfile = "output.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.getMass()
    print b.mass
    #toReaxLammps(b)

if __name__ == "__main__":
    pass
