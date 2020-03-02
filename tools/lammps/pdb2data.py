"""
@todo: assign atom types
"""
import os, sys
import numpy as np
import ase
import ase.io
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list

cutoff_table = {('H', 'H'):1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85, ('H', 'N'):1.3, ('N', 'N'):1.85, ('C', 'N'):1.85}
atom_types_table = {'C': 'opls_135', 'H': 'opls_140', 'N': 'opls_237'}
ERROR_ATP = 'Warning: User-defined atom types provided but the numbers do not match!\n'

def get_lists(atoms):

    # build neighbor list
    nl = []
    for i in range(len(atoms)):
        nl.append([])

    tokens_i, tokens_j = neighbor_list('ij', atoms, cutoff_table)
    for i in range(len(tokens_i)):
        nl[tokens_i[i]].append(tokens_j[i])

    # build bond list
    bond_list = []
    for i in range(len(nl)):
        if len(nl[i]) > 0:
            ai = i
            for j in range(len(nl[i])):
                aj = nl[i][j]
                if ai < aj:
                    bond_list.append([ai, aj])
    n_bond = len(bond_list)
    print(n_bond, "bond terms")

    angle_list = []
    # build angle list
    for i in range(len(nl)):
        if len(nl[i]) > 1:
            aj = i
            for j in range(len(nl[i])):
                ai = nl[i][j]
                for k in nl[i][j+1:]:
                    ak = k
                    angle_list.append([ai, aj, ak])
    n_angle = len(angle_list)
    print(n_angle, "angle terms")
    
    dihedral_list = []
    # build dihedral ist
    for i in bond_list:
        dj, dk = i[0], i[1]
        for j in nl[dj]:
            if j != dk:
                di = j
                for k in nl[dk]:
                    if k != dj:
                        dl = k
                        dihedral_list.append([di, dj, dk, dl])
    n_dihedral = len(dihedral_list)
    print(n_dihedral, "diheral terms")
    return nl, bond_list, angle_list, dihedral_list

def to_data(atoms, nl, bl, al, dl, bp, ap, dp):
    o = open('lammps.data', 'w')
    o.write('# lammps input file\n\n')
    o.write('%d atoms\n'%len(atoms))
    o.write('%d bonds\n'%len(bl))
    o.write('%d angles\n'%len(al))
    o.write('0 dihedrals\n')
    o.write('0 impropers\n')
    o.write('\n')

    o.write('1 atom types\n')
    o.write('1 bond types\n')
    o.write('1 angle types\n')
    o.write('0 dihedral types\n')
    o.write('0 improper types\n')
    o.write('\n')
    o.write('0.0 %12.7f xlo xhi\n'%atoms.cell[0][0])
    o.write('0.0 %12.7f ylo yhi\n'%atoms.cell[1][1])
    o.write('0.0 %12.7f zlo zhi\n'%atoms.cell[2][2])
    o.write('%12.7f %12.7f %12.7f xy xz yz\n'%(atoms.cell[1][0], atoms.cell[2][0], atoms.cell[2][1]))
    
    o.write('Masses\n\n')
    o.write('1 12.011\n')
    o.write('\n')

    o.write('Bond Coeffs\n\n')
    n = 1
    for i in bp:
        tokens = i.strip().split()
        o.write('%d %.6f %.4f\n'%(n, float(tokens[-2]), float(tokens[-1])))
        n += 1
    o.write('\n')

    o.write('Angle Coeffs\n\n')
    n = 1
    for i in ap:
        tokens = i.strip().split()
        o.write('%d %.6f %.4f\n'%(n, float(tokens[-2]), float(tokens[-1])))
        n += 1
    o.write('\n')
    #o.write('Dihedral Coeffs\n\n')

    n = 1
    o.write('Atoms\n\n')
    for i in atoms:
        o.write('%d 1 1 0.0'%n)
        for j in i.position:
            o.write('%10.5f'%j)
        o.write('\n')
        n += 1
    o.write('\n')

    n = 1
    o.write('Bonds\n\n')
    for i in bl:
        o.write('%d 1 %d %d\n'%(n, i[0]+1, i[1]+1))
        n +=1

    n = 1
    o.write('Angles\n\n')
    for i in al:
        o.write('%d 1 %d %d %d\n'%(n, i[0]+1, i[1]+1, i[2]+1))
        n +=1

    o.close()
"""
    n = 1
    o.write('Dihedrals\n\n')
    for i in dl:
        o.write('%d 1 %d %d %d %d\n'%(n, i[0]+1, i[1]+1, i[2]+1, i[3]+1))
        n +=1
"""


def read_ff():
    bp, ap, dp = [], [], []
    f = open('c.ff', 'r')
    for i in f:
        if 'BONDHARM:PARS' in i:
            bp.append(i)
        elif 'BENDAHARM:PARS' in i:
            ap.append(i)
        elif 'TORSION:PARS' in i:
            dp.append(i)
    f.close()
    return bp, ap, dp

if len(sys.argv) > 0:
    fname = sys.argv[1]
    atoms = ase.io.read(fname)
    nl, bl, al, dl = get_lists(atoms)
    bp, ap, dp = read_ff()
    to_data(atoms, nl, bl, al, dl, bp, ap, dp)

