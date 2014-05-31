"""
Generate a gaussian distribution of Al in CaO framework.

"""
import sys, os
import socket
import math
import ConfigParser, string
import copy
import numpy as np
import matplotlib.pyplot as plt
import random
import time

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"

sys.path.insert(0 , LIB)


from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo, toPdb
from cons import ATOMIC_MASS_CONSTANT as amc
from utilities import get_dist
from index import Group

def gaussian(x, delta, miu):
    f1 = -(x - miu)* (x-miu)
    f2 = f1/2/delta/delta
    f3 = 1/delta/math.sqrt(2*math.pi)*np.exp(f2)
    f4 = f3/np.sum(f3)
    f5 = f4*128
    f6 = np.around(f5/2) *2
    res = 128 - np.sum(f6)
    f6[15] += res/2
    f6[16] += res/2
    f7 = f6/2
    return f6, f7
    
def read_pdb():
    a = Pdb("out_sorted.pdb")
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    return b

def build_list(a):
    ca = []
    o = []

    n = 0
    m = 0
    ca.append([])
    o.append([])
    for i in a.atoms:
        if n > 0 and n%64 == 0:
            ca.append([])
            o.append([])
            m += 1
        if i.element == "Ca":
            ca[m].append(n) 
        elif i.element == "O":
            o[m].append(n)
        n += 1

    return ca, o

def replace(ca_list, o_list, al_natoms, o_natoms):
    al_list = []
    ignore_list = []
    for i in range(len(al_natoms)):
        tmp = random.sample(ca_list[i], al_natoms[i] + o_natoms[i])
        #print tmp, ca_list[i]
        if tmp:
            for j in range(len(tmp)):
                if j < al_natoms[i]:
                    al_list.append(tmp[j])
                else:
                    ignore_list.append(tmp[j])
        
    return al_list, ignore_list

def process(a, al_list, ignore_list):
    tmp = []
    for i in range(len(a.atoms)):
        if i in al_list:
            a.atoms[i].name = "Al"
            a.atoms[i].element = "Al"
            tmp.append(a.atoms[i])
        elif i in ignore_list:
            pass
        else:
            tmp.append(a.atoms[i])
    a.atoms = tmp

def main(counter):
    # Set up
    x = np.linspace(1,32,32)
    delta = 6.0
    miu = 17.0
    
    # get the distribution
    al_natoms, o_natoms = gaussian(x, delta, miu) 
    
    # read pdb
    a = read_pdb()
    
    # build the atoms list
    ca_list, o_list = build_list(a)
    
    # replace
    al_list, ignore_list = replace(ca_list, o_list, al_natoms, o_natoms)
    
    # do the replace
    process(a, al_list, ignore_list)
    
    # output the data
    toPdb(a, "case%02d.pdb"%counter)

    
if __name__ == "__main__":
    for i in range(20):
        main(i)
        time.sleep(3.0)
    
