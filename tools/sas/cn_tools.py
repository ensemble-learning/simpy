from operator import itemgetter

import sys, socket
import copy
import os
import math
import numpy as np

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
from output_conf import toPdb, toReaxLammps, toXyz

def read_data():
    datafile = "lammps.data"
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b

def surface_atoms_al(atoms, pbc):
    """
    Chemical Science.
    """
    sur = []
    r_coord_vector = []
    f = open("neighbourList.dat", "r")
    for i in f:
        tokens = i.strip().split(":")
        if len(tokens) > 0:
            xs = [0.0, 0.0, 0.0]
            rn = 0
            r2 = 0
            a0 = int(tokens[0])
            a1 = [int(j) for j in tokens[1].split()]
            x0 = atoms[a0-1][1:]
            #print a0, x0
            for k in a1:
                x2 = [0.0, 0.0, 0.0]
                x1 = atoms[k-1][1:]
                #print k, x1
                for l in range(3):
                    x2[l] = x1[l] - x0[l]
                for l in range(3):
                    if x2[l] > 0.5*pbc[l]: 
                        x2[l] = x2[l] - pbc[l]
                    elif x2[l] < -0.5*pbc[l]: 
                        x2[l] = x2[l] + pbc[l]
                    xs[l] += x2[l]
            for l in range(3):
                rn += xs[l]
                r2 += xs[l]*xs[l]
            r = math.sqrt(r2)
            r_coord_vector.append(r)

    return sur, r_coord_vector

def surface_atoms_cut_off(b):
    """
    """
    cn = np.loadtxt("cn.dat")
    o = open("sur.dat", "w")
    n_cut = 10
    sur = np.where(cn <= 10, 1, 0)
    
    b1 = copy.deepcopy(b)
    atoms = []
    for i in range(len(sur)):
        if sur[i] == 1:
            atoms.append(b.atoms[i])
            o.write("1\n")
        else:
            o.write("0\n")
    o.close()
    b1.atoms = atoms
    toPdb(b1, "surface_atoms.pdb")
    toXyz(b1, "surface_atoms.xyz")

def surface_atoms_read(b):
    sur = np.loadtxt("sur_sas.dat")
    b1 = copy.deepcopy(b)
    atoms = []
    for i in range(len(sur)):
        if int(sur[i]) == 1:
            atoms.append(b.atoms[i])
    b1.atoms = atoms
    toPdb(b1, "surface_atoms.pdb")
    toXyz(b1, "surface_atoms.xyz")

def main():
    b = read_data() 
    surface_atoms_method = 4
    if surface_atoms_method == 1:
        surface_atoms_al(atoms, pbc)
    elif surface_atoms_method == 2:
        sur = surface_atoms_cut_off(b)
    elif surface_atoms_method == 3:
        sur = surface_atoms_read(b)

if __name__ == "__main__":
    main()
