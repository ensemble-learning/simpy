import sys
import argparse

from mytype import System, Molecule, Atom
from poscar import Poscar
from output_conf import toXyz, toGeo, toPdb, toReaxLammps, toJdft
from cons import ATOMIC_MASS_CONSTANT

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="POSCAR", nargs="?", help="geo file name")
    args = parser.parse_args()
    
    poscar_file = args.fname
    a = Poscar(poscar_file)
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    vol = b.getVol()
    mass = b.getMass()
    density = mass/vol*ATOMIC_MASS_CONSTANT*1e27
    print density
    #print b.pbc
    
    toXyz(b, "contcar.xyz")
    toGeo(b, "geo")
    toPdb(b, "sim.pdb")
    toReaxLammps(b, "lammps.data")
    toJdft(b, "coords")


if __name__ == "__main__":
    main()
    
