#!/usr/bin/env python

from ase.io import read, write
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
from ase.build import molecule
import argparse

class Param():
    def __init__(self,):
        self.ncpu = '8'
        self.mem = '8gb'
        self.dft = 'b3lyp'
        self.basis = '6-311g(d,p)'
        self.solv = 'solvent=water, SMD'
        self.q = 0
        self.mult = 1

def gau_opt(p):
    calc_opt = Gaussian(label='opt', nprocshared=p.ncpu, mem=p.mem, save=None, method=p.dft, basis=p.basis, 
                                              charge=p.q, mult=p.mult)
    opt = GaussianOptimizer(atoms, calc_opt)
    opt.run()

    write('opt.xyz', atoms)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="g03.gjf", nargs="?", help="geometry file name")
    parser.add_argument("-q", type=int, help="charge")
    parser.add_argument("-nspin", type=int, help="Spin multiplicity")
    parser.add_argument("-ncpu", type=int, help="number of cpu for QM calculation")
    args = parser.parse_args()

    p = Param()
    infile = args.fname
    atoms = read(infile)

    if args.q:
        p.q = args.q
    if args.nspin:
        p.mult = args.nspin

    gau_opt(p)

