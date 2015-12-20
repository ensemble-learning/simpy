import sys
import argparse
from mytype import System, Molecule, Atom
from jaguar import Jaguar
from output_conf import toReaxLammps, toGeo, toPdb, toMsd, toXyz

def usage():
    print """python e_2_jag.py infile
    read the in file and output to geo and data files
    options:
    """

def convert(fname):
    a = Jaguar(fname)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    toXyz(b, "out.xyz")
    toPdb(b, "output.pdb", 1)
    toReaxLammps(b, "lammps.data")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="jag.in", nargs="?", help="geo file name")
    """
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    parser.add_argument("-sort", nargs=1, help="Sort the coordinations according to x, y or z")
    parser.add_argument("-element",action="store_true" , help="convert the atom name to element")
    """
    args = parser.parse_args()

    infile = args.fname
    convert(infile)

if __name__ == "__main__":
    main()
    


