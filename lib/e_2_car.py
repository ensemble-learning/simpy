""" manipulate .car file
"""

from mytype import System, Molecule, Atom
from car import Mdf, Car
from output_conf import toPdb

def sort_with_connect():
    """Sort the atoms: gather the same atoms, and atoms in same
    molecules according to connectivity (.mdf). This is very good 
    for generating gromacs input files, and organize the atoms into 
    a well defined order.
    """
    mdffile = "../testcode/car/LGPS_right_ms.mdf"
    carfile = "../testcode/car/LGPS_right_ms.car"
    mdf = Mdf(mdffile)
    car = Car(carfile)
    a = car.parser()
    Li = []
    Ge = []
    P = []
    atom_names = mdf.atom_names
    for i in range(len(mdf.atoms)):
        if mdf.atoms[i].element == "Li":
            Li.append(i)
        elif mdf.atoms[i].element == "Ge":
            Ge.append(i)
            for j in mdf.atoms[i].connect:
                Ge.append(atom_names.index(j))
        elif mdf.atoms[i].element == "P":
            P.append(i)
            for j in mdf.atoms[i].connect:
                P.append(atom_names.index(j))
    seq_new = Li + Ge + P
    a.sortNdx(seq_new)
    toPdb(a)

if __name__ == "__main__":
    sort_with_connect()
