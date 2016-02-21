""" Here we need to cut half of the nano wire to better view
"""
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
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/lib"
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from data import ReaxData
from output_conf import toPdb, toReaxLammps

def read_data():
    datafile = "lammps.data"
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b

def wrap(b):
    pbc = b.pbc
    for i in range(len(b.atoms)):
        for j in range(3):
            b.atoms[i].x[j] = b.atoms[i].x[j] - math.floor(b.atoms[i].x[j]/pbc[j])*pbc[j]

def get_center_along_z(b):
    dx = 2.77*2
    nl = int(b.pbc[2]/dx)
    dx = b.pbc[2]/nl
    
    center = []
    center_natoms = [0]*nl

    for i in range(nl):
        center.append([0.0, 0.0, 0.0])

    for i in range(len(b.atoms)):
        z = b.atoms[i].x[2]
        zn = int(z/dx)
        for j in range(3):
            center[zn][j] += b.atoms[i].x[j]
        center_natoms[zn] += 1
    for i in range(nl):
        if center_natoms[i] > 0:
            for j in range(3):
                center[i][j] = center[i][j]/center_natoms[i]
    return center, dx, nl

def filter_center(center, dx, nl, b):
    o = open("ndx_original.dat", "w")
    atoms = []
    for i in range(len(b.atoms)):
        x = b.atoms[i].x[0]
        z = b.atoms[i].x[2]
        an = b.atoms[i].an
        zn = int(z/dx)
        if x >= center[zn][0]:
            atoms.append(b.atoms[i])
            o.write("%d\n"%an)
    o.close()
    b.atoms = atoms
    toPdb(b)

        
def main():
    b = read_data()
    wrap(b)
    center, dx, nl, = get_center_along_z(b)
    filter_center(center, dx, nl, b)

main()
    
