""" read the pdb file and output to data (LAMMPS).
"""

import sys
import argparse
from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo, toPdb, toMsd, toXyz, toMusic, toPoscar, toJdft, toCfg, toTowheecoords

def usage():
    print("""python e_2_pdb [pbc|nopbc]
    read the pdb file and output to geo and data files
    options:
    pbc      crystal strucutre (default is supper.pdb)
    nobbc    gas phase cluster (default is output.pdb)
    """)

def test():
    testfile = "e_2_pdb.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.pbc = [20.80, 20.80, 20.80, 90.0, 90.0, 90.0]
    toReaxLammps(b)

def fortranOut(testfile = "output.pdb", args=''):
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    b.pbc = [50.00, 50.00, 50.00, 90.0, 90.0, 90.0]
    b.geotag = "XTLGRF 200"
    toReaxLammps(b, "lammps.data")
    if args:
        if args.element:
            toPdb(b, "out.pdb", 1)
    else:
        toPdb(b, "out.pdb")
    toMsd(b, "sim.msd")
    toGeo(b, "sim.geo")

def withPbc(testfile="supper.pdb", args=''):
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    b.toFrac()
    #b.translate(12.0, "z")
    toXyz(b, "out.xyz")
    toMusic(b, "out.music")
    if len(b.pbc) == 0:
        b.geotag = "BIOGRF 200"
    else:
        b.geotag = "XTLGRF 200"
    toReaxLammps(b, "lammps.data")
    toGeo(b, "sim.geo")
    toMsd(b, "sim.msd")
    toPoscar(b,)
    toCfg(b,)
    toJdft(b,)
    toTowheecoords(b)
    if args:
        if args.element:
            toPdb(b, "out.pdb", 1)
    else:
        toPdb(b, "out.pdb")

def sortXYZ(testfile="input.pdb", axis="z"):
    a = Pdb(testfile)
    b = a.parser()
    atoms = b.atoms
    b.sortXYZ(axis)
    toPdb(b, "out_sorted.pdb")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="input.pdb", nargs="?", help="geo file name")
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    parser.add_argument("-sort", nargs=1, help="Sort the coordinations according to x, y or z")
    parser.add_argument("-element",action="store_true" , help="convert the atom name to element")
    args = parser.parse_args()
    
    pdbfile = args.fname
    
    if args.c:
        if args.pbc:
            fortranOut(pdbfile, args)
        else:
            withPbc(pdbfile, args)

    if args.sort:
        sortXYZ(pdbfile, args.sort[0])
    
if __name__ == "__main__":
    main()
