""" 
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
    a = a_ref
    atoms = crystal([ele2, ele1], 
                [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)], 
                spacegroup=224, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C3"):
        os.mkdir("C3")
    os.chdir("C3")
    write_vasp("POSCAR", atoms, label="C3", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # C4 TiO2 (rutile) 
    a = a_ref
    c = a_ref * 0.6440
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.0, 0.0), (0.3057, 0.3057, 0.0)], 
                spacegroup=136, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C4"):
        os.mkdir("C4")
    os.chdir("C4")
    write_vasp("POSCAR", atoms, label="C4", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

    # C5 TiO2 (anatase) 
    a = a_ref
    c = a_ref * 2.9
    atoms = crystal([ele1, ele2], 
                [(0.0, 0.0, 0.0), (0.0, 0.0, 0.2081)], 
                spacegroup=141, cellpar=[a, a, a, 90, 90, 90])
    if not os.path.exists("C5"):
        os.mkdir("C5")
    os.chdir("C5")
    write_vasp("POSCAR", atoms, label="C5", direct=True, long_format=True, vasp5=True)
    os.chdir("..")

def main():
    generate_vaspinput()

if __name__ == "__main__":
    main()
