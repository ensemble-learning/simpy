"""
"""
import sys, os
import random
import ConfigParser, string
import socket

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from index import Group

def make_ndx(atp, lo, hi):
    assert lo < hi
    a = Pdb("model.pdb")
    b = a.parser()
    b.assignEleTypes()
    b.assignAtomTypes()
    n = 0
    counter = 0
    for i in b.atoms:
        if i.element == atp:
            z = i.x[2]
            if z > lo and z <= hi:
                print "%5d"%(n + 1),
                if (counter + 1) % 15 == 0:
                    print ''
                counter += 1
        n+= 1
    print ''

def main():
    print "[ Ca1 ]"
    make_ndx("Ca", 7.0, 12.0)
    print "[ Ca2 ]"
    make_ndx("Ca", 15.0, 19.0)
    print "[ Ca3 ]"
    make_ndx("Ca", 22.0, 27.0)
    print "[ Al1 ]"
    make_ndx("Al", 7.0, 12.0)
    print "[ Al2 ]"
    make_ndx("Al", 15.0, 19.0)
    print "[ Al3 ]"
    make_ndx("Al", 22.0, 27.0)

if __name__ == "__main__":

    main()


