"""
Generate inputs for scanning distance.
"""

import sys
import os
import socket

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

sys.path.insert(0 , LIB)

import argparse
import ConfigParser
from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toPoscar

"""
[SCAN]
points_values = 0.2850  0.3026  0.3201  0.3376  0.3551  0.3726  0.3901  0.4076  0.4251  0.4601  0.4952  0.5302  0.5652
ref_atom = 14
grp_atoms = 13 15 16
"""

class Scan():
    def __init__(self,):
        self.scan_points = []
        self.ref_atom = 0
        self.grp_atoms = []
        self.n_points = 0

def ConfigSectionMap(config, section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def read_conf(scan):

    config = ConfigParser.ConfigParser()
    config.read("scan.ini")

    tokens = ConfigSectionMap(config, "SCAN")['points_values']
    scan.scan_points = [float(i) for i in tokens.strip().split()]
    
    tokens = ConfigSectionMap(config, "SCAN")['ref_atom']
    scan.ref_atom = int(tokens)
    print scan.ref_atom
    
    tokens = ConfigSectionMap(config, "SCAN")['grp_atoms']
    scan.grp_atoms = [int(i) for i in tokens.strip().split()]
    print scan.grp_atoms


def gen_input(scan, poscar):
    n = len(scan.scan_points)
    n_ref = scan.ref_atom - 1
    for i in range(n):
        folder = "s%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        z0 = poscar.atoms[n_ref].xFrac[2]
        z1 = scan.scan_points[i]
        dz = z1 - z0
        # update the coords
        poscar.atoms[n_ref].xFrac[2] += dz 
        for j in range(len(scan.grp_atoms)):
            n_grp = scan.grp_atoms[j] - 1 
            poscar.atoms[n_grp].xFrac[2] += dz 
        # output the coords
        toPoscar(poscar)
        
        os.chdir("..")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="POSCAR", nargs="?", help="geo file name")
    args = parser.parse_args()
    
    scan = Scan()
    
    poscar_file = args.fname
    a = Poscar(poscar_file)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()

    # Read Config
    read_conf(scan)

    # Generate input
    gen_input(scan, b)

if __name__ == "__main__":
    main()


