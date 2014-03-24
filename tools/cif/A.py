""" Generate QM input (vasp) for A1(FCC), A2(BCC), A3(HCP), A4(diamond) and Ah (simple cubic).
Usage: python A.py element distance
Hint: distance reference to HCP structure, which is usually need scale 1.13 to BCC and 1.42 to FCC
"""

import os, shutil
import sys
from ase.lattice.spacegroup import crystal
from ase.io import read, write
from ase.io.vasp import write_vasp

def generate_vaspinput(element="Li", a_ref=3.05):

    # HCP
    a = a_ref
    c = 1.62 * a_ref
    atoms = crystal(element, [(1./3., 2./3., 3./4.)], spacegroup=194,
                cellpar=[a, a, c, 90, 90, 120])
    if not os.path.exists("HCP"):
        os.mkdir("HCP")
    os.chdir("HCP")
    write_vasp("POSCAR", atoms, label="HCP", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # FCC structure
    a = a_ref * 1.42 
    atoms = crystal(element, [(0,0,0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("FCC"):
        os.mkdir("FCC")
    os.chdir("FCC")
    write_vasp("POSCAR", atoms, label="FCC", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # BCC structure
    a = a_ref * 1.13
    atoms = crystal(element, [(0,0,0)], spacegroup=229, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("BCC"):
        os.mkdir("BCC")
    os.chdir("BCC")
    write_vasp("POSCAR", atoms, label="BCC", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # SCC structure
    a = a_ref * 0.90
    atoms = crystal(element, [(0,0,0)], spacegroup=221, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("SCC"):
        os.mkdir("SCC")
    os.chdir("SCC")
    write_vasp("POSCAR", atoms, label="SCC", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # Diamond
    a = a_ref * 1.93
    atoms = crystal(element, [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("DIAMOND"):
        os.mkdir("DIAMOND")
    os.chdir("DIAMOND")
    write_vasp("POSCAR", atoms, label="DIAMOND", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

def main():
    print __doc__
    if len(sys.argv) > 2:
        element = sys.argv[1]
        distance = float(sys.argv[2])
        generate_vaspinput(element, distance)

if __name__ == "__main__":
    main()
