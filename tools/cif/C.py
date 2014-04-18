""" 
Generate VASP input for C1, C2, C3, C4 and C5.
python C.py element1(Ti) element2(O2) a_ref.
We take C1 (CaF2) as reference. Usually C2 = C1, C3 = 0.9170* C1, C4 = 0.9278*C1 and C5 = 0.7798*C1.
"""

import os, shutil
import sys
from ase.lattice.spacegroup import crystal
from ase.io import read, write
from ase.io.vasp import write_vasp

def generate_vaspinput(ele1="Ti", ele2="O", a_ref=5.05):

    # C1 CaF2
    a = a_ref
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.5, 0.5), (0.25, 0.25, 0.25)], 
                spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C1"):
        os.mkdir("C1")
    os.chdir("C1")
    write_vasp("POSCAR", atoms, label="C1", direct=True, long_format=True, vasp5=True)
    os.chdir("..")
    
    # C2 FeS2
    a = a_ref
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.0, 0.0), (0.38504, 0.38504, 0.38504)], 
                spacegroup=205, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C2"):
        os.mkdir("C2")
    os.chdir("C2")
    write_vasp("POSCAR", atoms, label="C2", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # C3 Cu2O
    a = a_ref * 0.917
    atoms = crystal([ele2, ele1], 
                [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)], 
                spacegroup=224, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C3"):
        os.mkdir("C3")
    os.chdir("C3")
    write_vasp("POSCAR", atoms, label="C3", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # C4 TiO2 (rutile) 
    a = a_ref * 0.9278
    c = a * 0.6440
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.0, 0.0), (0.3057, 0.3057, 0.0)], 
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
    if not os.path.exists("C4"):
        os.mkdir("C4")
    os.chdir("C4")
    write_vasp("POSCAR", atoms, label="C4", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # C5 TiO2 (anatase) 
    a = a_ref * 0.7798
    c = a * 2.514
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.0, 0.0), (0.0, 0.0, 0.2081)], 
                spacegroup=141, cellpar=[a, a, c, 90, 90, 90])
    if not os.path.exists("C5"):
        os.mkdir("C5")
    os.chdir("C5")
    write_vasp("POSCAR", atoms, label="C5", direct=True, long_format=True, vasp5=True)
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
