""" parse the geo file with multi configuration into
seperated files
"""
import os
import argparse
from utilities import parseBlock
from mytype import System, Molecule, Atom
from geo import Geo
from output_conf import toGeo, toXyz, toReaxLammps, toPdb, toGjf

def convertors(b):
    toGeo(b, b.name+'.geo')
    toXyz(b, b.name+'.xyz')
    toXyz(b, b.name+'.gjf')
    toReaxLammps(b)
    toPdb(b)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="geo", nargs="?", help="geo file name")
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    parser.add_argument("-b", nargs=2, type=int, help="get the bond distance between a1, a2, a3")
    parser.add_argument("-a", nargs=3, type=int,help="get the angle of a1-a2-a3")
    parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")
    args = parser.parse_args()
    #print b.getBondDist(3,2)
    fname = args.fname

    assert os.path.exists(fname)
    a = Geo(fname)
    b = a.parser()
    b.assignAtomTypes()

    if args.c:
        print "converting %s to geo, xyz, gjf and lammps..."%fname
        convertors(b)

    if args.b:
        a1 = args.b[0] 
        a2 = args.b[1]
        val = b.getBondDist(a1, a2)
        print "Distance between %d and %d is %.3f."%(a1, a2, val)

    if args.a:
        a1 = args.a[0] 
        a2 = args.a[1]
        a3 = args.a[2]
        val = b.getAngle(a1, a2, a3)
        print "Angle of %d-%d-%d is %.3f."%(a1, a2, a3, val)

    if args.vol:
        vol = b.getVol()
        print "Volume is %.3f"%vol

if __name__ == "__main__":
    main()


