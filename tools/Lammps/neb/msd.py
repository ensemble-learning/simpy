"""
Calculate the displacement of atoms between t0 and
t1 state.

"""

import sys, math, os
import socket
import argparse

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

from dump import Dump
from output_conf import toXyz, toPdb

def gen_msd(log, dr):

    a0 = Dump("t0.lammpstrj")
    b0 = a0.parser()

    a1 = Dump("t1.lammpstrj")
    b1 = a1.parser()

    n0 = len(b0.atoms)
    n1 = len(b1.atoms)
    log.write

    log.write("t0 has %d atoms\n"%n0)
    log.write("t1 has %d atoms\n"%n1)

    if n0 != n1:
        sys.stderr.write("Atom numbers are different!")
        sys.stderr.flush()
        sys.exit()
    else:
        natom = n0

    neb_atoms = []
    o = open("detail.dat", "w")
    for i in range(natom):
        r2 = 0
        for j in range(3):
            d = (b0.atoms[i].x[j] - b1.atoms[i].x[j])
            r2 += d * d
        if r2 > dr*dr:
            neb_atoms.append(b1.atoms[i])
            o.write("%d\t%.4f\t%.4f\t%.4f\t%.4f\n"
                    %(i+1, math.sqrt(r2), b1.atoms[i].x[0],
                      b1.atoms[i].x[1], b1.atoms[i].x[2]))
    o.close()
    o = open("final.coords", "w")
    o.write("%d\n"%len(neb_atoms))
    for i in neb_atoms:
        o.write("%d\t%.4f\t%.4f\t%.4f\n"
                %(i.an, i.x[0], i.x[1], i.x[2]))
    o.close()

def main(args):

    log = open("gen_neb.log", "w")

    dr = 0.3
    if args.dr:
        dr = args.dr[0]
    log.write("Min displacement using is %.4f\n"%dr)

    gen_msd(log, dr)
    sys.stderr.close()
    log.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-dr", nargs=1, type=float, help="min displacement")
    args = parser.parse_args()
    main(args)
