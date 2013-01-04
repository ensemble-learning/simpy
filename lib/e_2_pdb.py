""" read the pdb file and output to data (LAMMPS).
"""
from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo

def test():
    testfile = "e_2_pdb.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.pbc = [20.80, 20.80, 20.80, 90.0, 90.0, 90.0]
    toReaxLammps(b)

def main():
    testfile = "supper.pdb"
    a = Pdb(testfile)
    b = a.parser()
    #b.assignAtomTypes()
    b.geotag = "XTLGRF 200"
    toReaxLammps(b)
    toGeo(b, "supper.geo")

if __name__ == "__main__":
    main()
