""" lammps force field file format parser
"""
import os
import math
import argparse
DEBUG = 0

class Ffield():
    """ lammps force field (ffield) file
    """
    def __init__(self, filename="ffield", lg=0):
        self.params = []
        self.pair = []
        self.bond = []
        self.angle = []
        self.torsion = []
        self.read(filename)
        """@xxxxx: read the input file """

    def read(self, filename):
        """ read the ffield file
        """
        f = open(filename, 'r')
        for i in f:
            if len(i.strip())> 0:
                if '#' in i:
                    tokens = i.strip().split('#')
                    data = tokens[0].strip().split()
                    comments = tokens[1]
                else:
                    data = i.strip().split()
                    comments = ''

                data.append(comments)
                if data[0] == 'pair_coeff':
                    self.pair.append(data)
                if data[0] == 'bond_coeff':
                    self.bond.append(data)
                if data[0] == 'angle_coeff':
                    self.angle.append(data)

        self.pair.sort(key=lambda x: x[1])
        self.bond.sort(key=lambda x: x[1])
        self.angle.sort(key=lambda x: x[1])

        f.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="ffield", nargs="?", help="force field file name")
    parser.add_argument("-D", action="store_true", help="Debug the code")
    parser.add_argument("-trans", action="store_true", help="transform ffield to more readable format")
    parser.add_argument("-params", action="store_true", help="generate params for training")
    parser.add_argument("-complete", action="store_true", help="complete the off table")
    parser.add_argument("-type", nargs=1, type=int, help="Force field type: 0 for vdw; 1 for lg_inner wall")
    parser.add_argument("-checkout", nargs='+', help="check out the force field")
    parser.add_argument("-check", action="store_true", help="check the force field")
    args = parser.parse_args()
    #print(b.getBondDist(3,2)
    
    fname = args.fname

    assert os.path.exists(fname)

    if args.D:
        DEBUG = 1

    if args.type:
        ntype = args.type[0]
    else:
        ntype = 0
        print("Note: Using default force field type")

    ff = Ffield(fname, ntype)
    if args.trans:
        ff.toEquation()
    if args.params:
        ff.toParams()
    if args.complete:
        ff.completeOff()
        toFfield(ff)
    if args.checkout:
        atoms = args.checkout
        ff2 = ff.checkout(atoms)
        toFfield(ff2)
    if args.check:
        ff.checkRedudant()

def test1():
    fname = "ffield"
    ntype = 0
    ff = Ffield(fname, ntype)
    ff2 = ff.checkout()
    toFfield(ff2)

def test():
    fname = "ffield"
    ntype = 0
    ff = Ffield(fname, ntype)
    ff.checkRedudant()
    
if __name__ == "__main__":
    #test()
    main()
