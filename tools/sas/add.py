""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import sys, socket
import copy
import os
import math
import numpy as np
import random

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
from pdb import Pdb
from output_conf import toPdb, toReaxLammps

DELTA_X = 2.77


def read_pdb(pdbfile):
    testfile = "out_sorted.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    return b

def read_cn(cnfile):
    cn = np.loadtxt(cnfile)
    return cn
    
def cal_centers(b):
    """
    Calculate the centers
    """
    dx = DELTA_X
    zl = b.pbc[2]
    n_slab = int(zl/dx) + 1
    dx = zl/n_slab
    
    n_atoms = [0]*n_slab
    centers = []
    for i in range(n_slab):
        centers.append([0.0, 0.0, 0.0])
    
    for i in b.atoms:
        x = i.x[0]
        y = i.x[1]
        z = i.x[2]
        # pbc
        z = z - int(z/zl)*zl
        n = int(z/dx)
        n_atoms[n] += 1
        centers[n][0] += x
        centers[n][1] += y
        centers[n][2] += z
    
    for i in range(n_slab):
        centers[i][0] = centers[i][0]/n_atoms[i]
        centers[i][1] = centers[i][1]/n_atoms[i]
        centers[i][2] = centers[i][2]/n_atoms[i]
    return centers

def add_oh(b, centers, cn):
    dx = DELTA_X
    zl = b.pbc[2]
    b1 = copy.deepcopy(b)
    n_cut = 10
    for i in range(len(cn)):
        print "processing %d"%(i*1.0/8183*100)
        token = random.random()
        n = int(cn[i])
        if n < n_cut and token > 0.9:
            a1 = Atom()
            a1.name = "O"
            a1.element = "O"
            x0 = b.atoms[i].x[0]
            y0 = b.atoms[i].x[1]
            z0 = b.atoms[i].x[2]
            z0 = z0 - int(z0/zl)*zl
            n_sl = int(z0/dx)
            x1 = centers[n_sl][0]
            y1 = centers[n_sl][1]
            z1 = centers[n_sl][2]
            r2 = (x0-x1)*(x0-x1)
            r2 += (y0-y1)*(y0-y1)
            r2 += (z0-z1)*(z0-z1)
            r1 = math.sqrt(r2)
            r = r1 + 2.0
            s = r/r1
            x = s*(x0-x1) + x1 + 1.0
            y = s*(y0-y1) + y1 + 1.0
            z = s*(z0-z1) + z1
            a1.x[0] = x 
            a1.x[1] = y
            a1.x[2] = z
            b1.atoms.append(a1)

            a2 = Atom()
            a2.name = "H"
            a2.element = "H"
            r = r1 + 3.0
            s = r/r1
            x = s*(x0-x1) + x1 + 1.0
            y = s*(y0-y1) + y1 + 1.0
            z = s*(z0-z1) + z1
            a2.x[0] = x 
            a2.x[1] = y
            a2.x[2] = z
            b1.atoms.append(a2)

    toPdb(b1, "add_oh.pdb")

def add_o(b, centers, cn):
    dx = DELTA_X
    zl = b.pbc[2]
    b1 = copy.deepcopy(b)
    n_cut = 10
    for i in range(len(cn)):
        print "processing %d"%(i*1.0/8183*100)
        n = int(cn[i])
        token = random.random()
        if n < n_cut and token > 0.9:
            a = Atom()
            a.name = "O"
            a.element = "O"
            x0 = b.atoms[i].x[0]
            y0 = b.atoms[i].x[1]
            z0 = b.atoms[i].x[2]
            z0 = z0 - int(z0/zl)*zl
            n_sl = int(z0/dx)
            x1 = centers[n_sl][0]
            y1 = centers[n_sl][1]
            z1 = centers[n_sl][2]
            r2 = (x0-x1)*(x0-x1)
            r2 += (y0-y1)*(y0-y1)
            r2 += (z0-z1)*(z0-z1)
            r1 = math.sqrt(r2)
            r = r1 + 2.0
            s = r/r1
            x = s*(x0-x1) + x1 + 1.0
            y = s*(y0-y1) + y1 + 1.0
            z = s*(z0-z1) + z1
            a.x[0] = x 
            a.x[1] = y
            a.x[2] = z
            b1.atoms.append(a)
    toPdb(b1, "add_o.pdb")

def add_h(b, centers, cn):
    dx = DELTA_X
    zl = b.pbc[2]
    b1 = copy.deepcopy(b)
    n_cut = 10
    for i in range(len(cn)):
        n = int(cn[i])
        print "processing %d"%(n*1.0/8183*100)
        if n < n_cut:
            a = Atom()
            a.name = "H"
            a.element = "H"
            x0 = b.atoms[i].x[0]
            y0 = b.atoms[i].x[1]
            z0 = b.atoms[i].x[2]
            z0 = z0 - int(z0/zl)*zl
            n_sl = int(z0/dx)
            x1 = centers[n_sl][0]
            y1 = centers[n_sl][1]
            z1 = centers[n_sl][2]
            r2 = (x0-x1)*(x0-x1)
            r2 += (y0-y1)*(y0-y1)
            r2 += (z0-z1)*(z0-z1)
            r1 = math.sqrt(r2)
            r = r1 + 1.5
            s = r/r1
            x = s*(x0-x1) + x1
            y = s*(y0-y1) + y1
            z = s*(z0-z1) + z1
            a.x[0] = x 
            a.x[1] = y
            a.x[2] = z
            b1.atoms.append(a)
    toPdb(b1, "add_h.pdb")
    #toPdb(b1, "add_h_%05d.pdb"%i)
    
def main():
    pdbfile = "out_sorted.pdb"
    cnfile = "cn.aux"
    b = read_pdb(pdbfile)
    centers = cal_centers(b)
    cn = read_cn(cnfile)
    #add_h(b, centers, cn)
    #add_o(b, centers, cn)
    add_oh(b, centers, cn)
    

main()
