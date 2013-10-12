""" manipulate .car file
@log: 
2013_09_29: add a function __get_index() to get the atom index
"""

import argparse
from mytype import System, Molecule, Atom
from car import Mdf, Car
from output_conf import toPdb

def __get_index(atom_names, atom):
    """ Parse the atom index according to the atom name.
    @note: still not sure how ms works in assigning the 
    atom name.
    """
    n = 0
    for i in atom_names:
        if atom in i:
            break
        n += 1
    return n

def sort_with_connect(fname):
    """Sort the atoms: gather the same atoms, and atoms in same
    molecules according to connectivity (.mdf). This is very good 
    for generating gromacs input files, and organize the atoms into 
    a well defined order.
    """
    mdffile = "%s.mdf"%fname
    carfile = "%s.car"%fname
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
                P.append(__get_index(atom_names, j))
    seq_new = Li + Ge + P
    a.sortNdx(seq_new)
    toPdb(a)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="input", nargs="?", help="car file name (without .car)")
    parser.add_argument("-sort", action="store_true", help="sort the coordination according to connectivities")
    args = parser.parse_args()
    fname = args.fname
    if args.sort:
        sort_with_connect(fname)
