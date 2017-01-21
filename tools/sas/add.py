""" add adsorbed atoms on catalysis
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
elif "tao-ThinkCentre-M79" in socket.gethostname():
    LIB = "/home/tao/Soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from data import ReaxData
from output_conf import toPdb, toReaxLammps

DELTA_X = 2.77

def read_pdb(pdbfile):
    """read pdb file
    """

    a = Pdb(pdbfile)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    return b

def read_data(datafile):
    """read data file
    """
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    return b

def read_cn(cnfile):
    """read cn file
    """
    cn = np.loadtxt(cnfile)
    return cn

def read_build_log():
    """read log file 
    """
    n_ignore = 0
    ndx = []
    f =open("build.log", "r")
    for i in f:
        if i.strip().startswith("Fixed atoms is"):
            tokens = i.strip().split()
            n_ignore = int(tokens[-1].strip("."))
        elif i.strip().startswith("Index is as following:"):
            break
    for i in f:
        if i.strip().startswith("Output to pdb."):
            break
        else:
            tokens = i.strip().split()
            if len(tokens) > 0:
                ndx += [int(j) for j in tokens]
    f.close()
    print ndx
    return n_ignore, ndx
    
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

def add_oh(b, centers, cn, n_ignore, ndx):
    """add oh
    """
    dx = DELTA_X
    zl = b.pbc[2]
    b1 = copy.deepcopy(b)
    b2 = copy.deepcopy(b)
    
    a = Atom()
    a.name = "O"
    a.element = "O"
    b2.atoms.append(a)
    a = Atom()
    a.name = "H"
    a.element = "H"
    b2.atoms.append(a)
    
    n_cut = 10
    r_pt_oh = 2.0
    r_oh = 1.0
    for i in range(len(cn)):
        #print "processing %d"%(i*1.0/len(cn)*100)
        token = random.random()
        n = int(cn[i])
        if n < n_cut and token >= 0.0:
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
            r = r1 + r_pt_oh
            s = r/r1
            x = s*(x0-x1) + x1 
            y = s*(y0-y1) + y1 
            z = s*(z0-z1) + z1
            a1.x[0] = x 
            a1.x[1] = y
            a1.x[2] = z
            b1.atoms.append(a1)
            b2.atoms[-2].x[0] = x
            b2.atoms[-2].x[1] = y
            b2.atoms[-2].x[2] = z

            a2 = Atom()
            a2.name = "H"
            a2.element = "H"
            r = r1 + r_pt_oh + r_oh
            s = r/r1
            x = s*(x0-x1) + x1 
            y = s*(y0-y1) + y1 
            z = s*(z0-z1) + z1
            a2.x[0] = x 
            a2.x[1] = y
            a2.x[2] = z
            b1.atoms.append(a2)
            b2.atoms[-1].x[0] = x
            b2.atoms[-1].x[1] = y
            b2.atoms[-1].x[2] = z
            if i >= n_ignore - 1:
                toPdb(b2, "%05d_%05d.pdb"%(i, ndx[i]))

    toPdb(b1, "add_oh.pdb")

def add_o(b, centers, cn):
    """add o
    """
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
    """add h
    """
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
    """add adsorbed atoms on catalysis
    """
    pdb = 0
    if pdb:
        pdbfile = "out_sorted.pdb"
        b = read_pdb(pdbfile)
    data = 1
    if data:
        datafile = "lammps.data"
        b = read_data(datafile)
        
    centers = cal_centers(b)

    cnfile = "cn.aux"
    if os.path.exists("build.log"):
        n_ignore, ndx = read_build_log()
    cn = read_cn(cnfile)
    #add_h(b, centers, cn)
    #add_o(b, centers, cn)
    add_oh(b, centers, cn, n_ignore, ndx)
    

if __name__ == "__main__":
    main()
