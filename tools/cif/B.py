""" 
Generate B1, B2, B3 and B4 structures.
python B.py element1 element2 r_ref
Here we use B1 strucutre as r_ref. Usually, B2 = 0.6975 * B1, B3 = 1.086*B1 and B4 = 0.8272*B1.
"""

import os, shutil
import sys
from ase.lattice.spacegroup import crystal
from ase.io import read, write
from ase.io.vasp import write_vasp

def generate_vaspinput(ele1="Ca", ele2="O", a_ref=3.05):

    # B1 NaCl
    a = a_ref
    atoms = crystal([ele1, ele2], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("B1"):
        os.mkdir("B1")
    os.chdir("B1")
    write_vasp("POSCAR", atoms, label="B1", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # B2 CsCl
    a = a_ref * 0.6975
    atoms = crystal([ele1, ele2], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=221,
               cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("B2"):
        os.mkdir("B2")
    os.chdir("B2")
    write_vasp("POSCAR", atoms, label="B2", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # B3 ZnS Sphalerite (F43m)
    a = a_ref * 1.086
    atoms = crystal([ele1, ele2], [(0, 0, 0), (0.25, 0.25, 0.25)], spacegroup=216,
               cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("B3"):
        os.mkdir("B3")
    os.chdir("B3")
    write_vasp("POSCAR", atoms, label="B3", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # B4 ZnS Wurtzite (p6mc)
    a = a_ref * 0.8272
    c = a * 1.635
    atoms = crystal([ele1, ele2], [(0.3333, 0.6667, 0), (0.3333, 0.6667, 0.3748)], spacegroup=186,
               cellpar=[a, a, c, 90, 90, 120])
    if not os.path.exists("B4"):
        os.mkdir("B4")
    os.chdir("B4")
    write_vasp("POSCAR", atoms, label="B4", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

def main():
    
    if len(sys.argv) >= 4:
        ele1 = sys.argv[1]
        ele2 = sys.argv[2]
        ref = float(sys.argv[3])
        generate_vaspinput(ele1, ele2, ref)
    else:
        print __doc__

if __name__ == "__main__":
    main()
