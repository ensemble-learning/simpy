"""
Generate vasp isotropic scan
"""

import sys, os, copy
import socket
import shutil
import numpy as np
import argparse

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
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    LIB = "/home/tao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from poscar import Poscar
from utilities import lattice2v, v2lattice
from output_conf import toPoscar, toXyz

def gen_scan(s, args):

    start = 0.7
    end = 1.3
    start = np.power(start, 1/3.0)
    end = np.power(end, 1/3.0)
    n = 11
    x = np.linspace(start, end, n)

    for i in range(n):
        folder = "scan_%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        shutil.copy("../INCAR", ".")
        shutil.copy("../POTCAR", ".")
        shutil.copy("../KPOINTS", ".")
        if os.path.exists('../pbs'):
            shutil.copy("../pbs", ".")
        s_new = copy.copy(s)
        xx, xy, xz, yy, yz, zz = lattice2v(s.pbc)
        a = np.array([xx, 0.0, 0.0])
        b = np.array([xy, yy, 0.0])
        c = np.array([xz, yz, zz])
        if args.xyz:
            tmp = args.xyz[0] 
            if tmp == "xy":
                a = a * x[i]
                b = b * x[i]
            elif tmp == "z":
                c = c * x[i]
        else:
            a = a * x[i]
            b = b * x[i]
            c = c * x[i]
        pbc = v2lattice(a, b, c)
        s_new.pbc = pbc
        toPoscar(s_new, "POSCAR")
        os.chdir("..")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="POSCAR", nargs="?", help="geo file name")
    parser.add_argument("-xyz", nargs=1, help="xyz, xy, z")
    parser.add_argument("-scale", nargs=2, type=float, help="scale factor")
    args = parser.parse_args()
    #gen_scan(args)
    
    poscar_file = "POSCAR"
    a = Poscar(poscar_file)
    poscar_file = args.fname
    a = Poscar(poscar_file)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    gen_scan(b, args)

if __name__ == "__main__":
    main()
