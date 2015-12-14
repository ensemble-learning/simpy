""" read the geo file and output to data (LAMMPS), geo and xyz file.
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
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from data import ReaxData
from output_conf import toPdb, toReaxLammps

def read_data(datafile):
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b

def read_build_log(logfile):
    z0 = 0.0
    ndx = []
    f = open(logfile, "r")
    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("Starting from"):
            z0 = float(tokens[2])
        if i.strip().startswith("Index is as following:"):
            break
    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("Output to pdb."):
            break
        if len(tokens) > 0:
            ndx += [int(j) for j in tokens]
    return z0, ndx

def update_coords(d1, d2, z0, ndx):
    for i in range(len(ndx)):
        id = ndx[i] - 1
        d1.atoms[id].x[0] = d2.atoms[i].x[0]
        d1.atoms[id].x[1] = d2.atoms[i].x[1]
        d1.atoms[id].x[2] = d2.atoms[i].x[2] + z0
    
def main(log):
    d1_file = "lammps.data"
    d2_file = "./slice_00/conf2/lammps.data"
    d1 = read_data(d1_file)
    d2 = read_data(d2_file)
    logfile = "./slice_00/conf2/build.log"
    z0, ndx = read_build_log(logfile)
    toPdb(d1, "t0.pdb")
    update_coords(d1, d2, z0, ndx)
    toPdb(d1, "t1.pdb")
    toReaxLammps(d1, "lammps.data.1")

if __name__ == "__main__":
    log = open("gether.log", "w")
    main(log)
    log.close()
