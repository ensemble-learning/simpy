""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import sys, socket
import copy
import os
import math
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
from data import ReaxData
from output_conf import toPdb, toReaxLammps

def read_data():
    datafile = "lammps.data"
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b

data = read_data()

stotal = 0.0

pbc = data.pbc
XL = pbc[0]
YL = pbc[1]
ZL = pbc[2]
n_sample = 50
atomsigma = 4.2
pi = math.pi

atoms = data.atoms

for i in range(len(atoms)):
    print "Processing %d"%i
    ncount = 0
    for j in range(n_sample):
        phi = random.random()*2.0*pi
        costheta = 1 - random.random() * 2.0
        theta = math.acos(costheta)
        xpoint = math.sin(theta)*math.cos(theta)
        ypoint = math.sin(theta)*math.sin(theta)
        zpoint = costheta

        xpoint = xpoint*atomsigma/2.0
        ypoint = ypoint*atomsigma/2.0
        zpoint = zpoint*atomsigma/2.0
        
        xpoint = xpoint + atoms[i].x[0]
        ypoint = ypoint + atoms[i].x[1]
        zpoint = zpoint + atoms[i].x[2]

        if xpoint < 0.0:
            xpoint = xpoint + XL
        if xpoint >= XL:
            xpoint = xpoint - XL
        if ypoint < 0.0:
            ypoint = ypoint + YL
        if ypoint >= YL:
            ypoint = ypoint - YL
        if zpoint < 0.0:
            zpoint = zpoint + ZL
        if zpoint >= ZL:
            zpoint = zpoint - ZL

        deny = 0
        for k in range(len(atoms)):
            if k == i:
                pass
            else:
                dx = xpoint - atoms[k].x[0]
                dx = dx - XL*int(2.0*dx/XL)
                dy = ypoint - atoms[k].x[1]
                dy = dy - YL*int(2.0*dy/YL)
                dz = zpoint - atoms[k].x[2]
                dz = dz - ZL*int(2.0*dz/ZL)
                dist2 = dx*dx + dy*dy + dz*dz
                
                if math.sqrt(dist2) < 0.999*atomsigma/2.0:
                    deny = True
                    break
        if deny:
            pass
        else:
            ncount += 1
    sfrac = ncount*1.0/n_sample
    sjreal = pi*atomsigma*atomsigma*sfrac
    stotal += sjreal

print stotal
        
