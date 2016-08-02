""" read the geo file and output to data (LAMMPS), geo and xyz file.  """
import numpy as np
import copy
import argparse
from mytype import System, Molecule, Atom
from data import ReaxData, FullData
from output_conf import toPdb, toFullLammps

def extract_mols(b):
    mols = []
    nmols = b.atoms[-1].resn
    for i in range(nmols):
        mol = Molecule()
        mols.append(mol)
    
    for i in b.atoms:
        n = i.resn - 1
        mols[n].atoms.append(i)
        
    atoms_start = []
    atoms_end = []
    for i in mols:
        i.a0 = i.atoms[0].an
        i.a1 = i.atoms[-1].an
        i.n_atoms = i.a1 - i.a0 + 1
        atoms_start.append(i.a0)
        atoms_end.append(i.a1)
    
    for i in mols:
        if len(i.bonds) > 0:
            for j in i.bonds:
                j[0] = j[0] - i.b0 + 1
                j[2] = j[2] - i.a0 + 1
                j[3] = j[3] - i.a0 + 1
    
    for i in b.bonds:
        blist = [i[2], i[3]]
        blist.sort()
        for n in range(len(atoms_end)):
            if atoms_end[n] >= blist[1]:
                break
        mols[n].bonds.append(i)
    
    bonds_start = []
    bonds_end = []
    for i in mols:
        if len(i.bonds) > 0:
            i.b0 = i.bonds[0][0]
            i.b1 = i.bonds[-1][0]
            i.n_bonds = i.b1 - i.b0 + 1
            bonds_start.append(i.b0)
            bonds_end.append(i.b1)
    
    for i in mols:
        if len(i.bonds) > 0:
            for j in i.bonds:
                j[0] = j[0] - i.b0 + 1
                j[2] = j[2] - i.a0 + 1
                j[3] = j[3] - i.a0 + 1
    
    for i in b.angles:
        alist = [i[2], i[3], i[4]]
        alist.sort()
        for n in range(len(atoms_end)):
            if atoms_end[n] >= alist[2]:
                break
        mols[n].angles.append(i)
    
    angles_start = []
    angles_end = []
    for i in mols:
        if len(i.angles) > 0:
            i.ang0 = i.angles[0][0]
            i.ang1 = i.angles[-1][0]
            i.n_angles = i.ang1 - i.ang0 + 1
            angles_start.append(i.b0)
            angles_end.append(i.b1)
    
    for i in mols:
        if len(i.angles) > 0:
            for j in i.angles:
                j[0] = j[0] - i.ang0 + 1
                j[2] = j[2] - i.a0 + 1
                j[3] = j[3] - i.a0 + 1
                j[4] = j[4] - i.a0 + 1
    
    for i in b.dihedrals:
        alist = [i[2], i[3], i[4], i[5]]
        alist.sort()
        for n in range(len(atoms_end)):
            if atoms_end[n] >= alist[2]:
                break
        mols[n].dihedrals.append(i)
    
    dihedrals_start = []
    dihedrals_end = []
    for i in mols:
        if len(i.dihedrals) > 0:
            i.d0 = i.dihedrals[0][0]
            i.d1 = i.dihedrals[-1][0]
            i.n_dihedrals = i.d1 - i.d0 + 1
            dihedrals_start.append(i.d0)
            dihedrals_end.append(i.d1)
    
    for i in mols:
        if len(i.dihedrals) > 0:
            for j in i.dihedrals:
                j[0] = j[0] - i.d0 + 1
                j[2] = j[2] - i.a0 + 1
                j[3] = j[3] - i.a0 + 1
                j[4] = j[4] - i.a0 + 1
                j[5] = j[5] - i.a0 + 1
    return mols

def grow_system(b, mols, nmols):
    c = copy.copy(b)
    c.atoms = []
    c.bonds = []
    c.angles = []
    c.dihedrals = []

    n1, n2, n3, n4 = [0, 0, 0, 0]
    for i in range(len(nmols)):
        n1 += nmols[i]*mols[i].n_atoms 
        n2 += nmols[i]*mols[i].n_bonds
        n3 += nmols[i]*mols[i].n_angles
        n4 += nmols[i]*mols[i].n_dihedrals
    c.n_atoms = n1
    c.n_bonds = n2
    c.n_angles = n3
    c.n_dihedrals = n4

    n1, n2, n3, n4, n5 = [0, 0, 0, 0, 0]
    for i in range(len(nmols)): 
        n_repeat = nmols[i]
        for j in range(n_repeat):
            for k in mols[i].atoms:
                kn = copy.copy(k)
                kn.resn = n3 + 1
                c.atoms.append(kn)
                n2 += 1
            for k in mols[i].bonds:
                kn = copy.copy(k)
                kn[2] = k[2] + n1
                kn[3] = k[3] + n1
                c.bonds.append(kn)
            for k in mols[i].angles:
                kn = copy.copy(k)
                kn[2] = k[2] + n1
                kn[3] = k[3] + n1
                kn[4] = k[4] + n1
                c.angles.append(kn)
            for k in mols[i].dihedrals:
                kn = copy.copy(k)
                kn[2] = k[2] + n1
                kn[3] = k[3] + n1
                kn[4] = k[4] + n1
                kn[5] = k[5] + n1
                c.dihedrals.append(kn)
            n3 += 1
            n1 = n2
    return c

def fconv(fname):
    a = FullData(fname)
    b = a.parser()
    b.assignAtomTypes()
    toPdb(b)

def gen(fname, infile):
    f = open("inp", "r")
    for i in f:
        tokens = i.strip().split()
        nmols = [int(j) for j in tokens]
    print nmols
    a = FullData(fname)
    b = a.parser()
    b.assignAtomTypes()
    mols = extract_mols(b)
    c = grow_system(b, mols, nmols)
    c.assignAtomTypes()
    toFullLammps(c)

def update(fname):
    coords = np.loadtxt("coords")
    a = FullData(fname)
    b = a.parser()
    b.assignAtomTypes()
    for i in range(len(b.atoms)):
        b.atoms[i].x = coords[i]
    toFullLammps(b, "new.data")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="lammps.data", nargs="?", help="geo file name")
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    parser.add_argument("-gen", action="store_true", help="generate lammps input from template")
    parser.add_argument("-update", action="store_true", help="generate lammps input from template")
    parser.add_argument("-s", nargs=1, help="configuration file")
    args = parser.parse_args()

    datafile = args.fname
    if args.c:
        fconv(datafile)
    
    if args.gen:
        if args.s:
            infile = args.s[0]
        else:
            infile = "inp"
        gen(datafile, infile)

    if args.update:
        update(datafile)

if __name__ == "__main__":
    main()
