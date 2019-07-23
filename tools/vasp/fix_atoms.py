#!/usr/bin/env python

import sys, shutil
import argparse
import os, copy
import socket

LIB = ''

print(socket.gethostname())

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
elif "comet" in socket.gethostname():
    LIB = "/home/tcheng/soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "zwicky":
    LIB = "/home/tcheng/Soft/simpy/lib"
elif "onyx" in socket.gethostname():
    LIB = "/p/home/taocheng/src/simpy/lib"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    LIB = "/home/tao/soft/simpy/lib"
elif socket.gethostname() == "mu05":
    LIB = "/home/chengtao/soft/simpy/lib"
elif "stampede2" in socket.gethostname():
    LIB = "/home1/04076/tg833760/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toPoscar
from cons import ATOMIC_MASS_CONSTANT

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="POSCAR", nargs="?", help="geo file name")
    parser.add_argument("-n", type=int, help="number of fixed atoms")
    args = parser.parse_args()
    if args.n:
        n_fixed = int(args.n)
    else:
        n_fixed = int(input("number of fixed atoms: "))

    poscar_file = args.fname
    a = Poscar(poscar_file)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    vol = b.getVol()
    mass = b.getMass()
    density = mass/vol*ATOMIC_MASS_CONSTANT*1e27
    
    shutil.copy(poscar_file, poscar_file+".1")
    
    coords = []
    for i in range(len(b.atoms)):
        coords.append([])

    for i in range(len(b.atoms)):
        coords[i].append(i)
        for j in b.atoms[i].x:
            coords[i].append(j)
    coords.sort(key=lambda x: x[3])

    for i in range(len(coords)):
        if i < n_fixed:
            coords[i].append("F")
        else:
            coords[i].append("T")
    coords.sort(key=lambda x: x[0])

    for i in range(len(b.atoms)):
        b.atoms[i].xr = [0,0,0]
        if coords[i][-1] == "F":
            b.atoms[i].xr = [1,1,1]
    toPoscar(b, "POSCAR_new")
    

if __name__ == "__main__":
    main()
    
