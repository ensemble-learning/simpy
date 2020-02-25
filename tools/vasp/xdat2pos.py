#!/usr/bin/env python

import sys
import os
import socket
import copy

print(socket.gethostname())

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
elif socket.gethostname() == "zwicky":
    LIB = "/home/tcheng/Soft/simpy/lib"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/lib/"
elif socket.gethostname() == "fermion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif "onyx" in socket.gethostname():
    LIB = "/p/home/taocheng/src/simpy/lib"
elif "armstrong" in socket.gethostname():
    LIB = "/p/home/taocheng/src/simpy/lib"
elif "stampede2" in socket.gethostname():
    LIB = "/home1/04076/tg833760/soft/simpy/lib"
elif "login" in socket.gethostname():
    LIB = "/central/home/tcheng/soft/simpy/lib"
elif socket.gethostname() == "mu01":
    LIB = "/opt/sourcecoude/simpy/lib"
elif "tao-Precision-Tower-3420-ubuntu" in socket.gethostname():
    LIB = "/home/tao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toPoscar

def parse_XDATCAR(t0=0):
    head = []
    tag = 1
    f = open("XDATCAR", "r")
    for i in f:
        if i.strip().startswith("Direct configuration"):
            tokens = i.strip().split("=")
            step = int(tokens[1]) + t0
            break
        else:
            head.append(i)
    
    conf = []
    steps = []
    coords = []
    while(tag):
        tag = 0
        for i in f:
            tag = 1
            if i.strip().startswith("Direct configuration"):
                conf.append(coords)
                steps.append(step)
                coords = []
                tokens = i.strip().split("=")
                step = int(tokens[1]) + t0
                break
            else:
                coords.append(i)
    conf.append(coords)
    steps.append(step)
    
    if not os.path.exists("trj"):
        os.mkdir("trj")
    os.chdir("trj")
    
    for i in range(len(conf)):
        o = open("coord_%06d"%steps[i], "w")
        for j in conf[i]:
            o.write(j)
        o.close()
    os.chdir("..")
    return conf, steps

def toPOSCAR(t0, conf, steps):

    if not os.path.exists("poscars"):
        os.mkdir("poscars")
    os.chdir("poscars")

    for i in range(len(conf)):
        fname = "POSCAR_%06d"%steps[i]
        t0.name = "step = %d"%steps[i]
        for j in range(len(conf[i])):
            tokens = conf[i][j].strip().split()
            t0.atoms[j].xFrac[0] = float(tokens[0])
            t0.atoms[j].xFrac[1] = float(tokens[1])
            t0.atoms[j].xFrac[2] = float(tokens[2])
        toPoscar(t0, fname)

    os.chdir("..")

def main():
    t0 = 0
    if len(sys.argv) > 1:
        t0 = int(sys.argv[1])
    a = Poscar("POSCAR") 
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    conf, steps = parse_XDATCAR(t0)
    toPOSCAR(b, conf, steps)

if __name__ == "__main__":
    main()
        
