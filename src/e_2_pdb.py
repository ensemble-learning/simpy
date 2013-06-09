
""" read the pdb file and output to data (LAMMPS).
"""

import sys
sys.path.append("/home/tao/Nutstore/code/simupy/lib")

from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toReaxLammps, toGeo

def usage():
    print """python e_2_pdb [pbc|nopbc]
    read the pdb file and output to geo and data files
    options:
    pbc      crystal strucutre (default is supper.pdb)
    nobbc    gas phase cluster (default is output.pdb)
    """

def test():
    testfile = "e_2_pdb.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.pbc = [20.80, 20.80, 20.80, 90.0, 90.0, 90.0]
    toReaxLammps(b)

def fortranOut():
    testfile = "output.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.pbc = [50.00, 50.00, 50.00, 90.0, 90.0, 90.0]
    b.geotag = "XTLGRF 200"
    toReaxLammps(b, "lammps.data")
    toGeo(b, "sim.geo")

def withPbc():
    testfile = "supper.pdb"
    a = Pdb(testfile)
    b = a.parser()
    b.assignAtomTypes()
    b.geotag = "BIOGRF 200"
    toReaxLammps(b, "lammps.data")
    toGeo(b, "sim.geo")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "using default nopbc"
        fortranOut()
    else:
        opt = sys.argv[1]
        if opt == "pbc":
            withPbc()
        elif opt == "nopbc":
            fortranOut()
        else:
            usage()

