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
    testfile = "test.pdb"
    a = Pdb(testfile)
    b = a.parser()
    atoms = b.atoms
    print atoms[0].x
    b.sortZ()
    atoms = b.atoms
    print atoms[0].x

if __name__ == "__main__":
    test()
