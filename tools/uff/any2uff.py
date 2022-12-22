import os
from ase.io import read, write
import json

uff_atp = {"H": "H_", "He": "He4+4", "Li": "Li", "Be": "Be3+2", "B": "B_2", "C": "C_3", "N": "N_3", "O": "O_3", "F": "F_", "Ne": "Ne4+4", "Na": "Na", "Mg": "Mg3+2", "Al": "Al3", "Si": "Si3", "P": "P_3+3", "S": "S_3+2", "Cl": "Cl", "Ar": "Ar4+4", "K": "K_", "Ca": "Ca6+2", "Sc": "Sc3+3", "Ti": "Ti6+4", "V": "V_3+5", "Cr": "Cr6+3", "Mn": "Mn6+2", "Fe": "Fe6+2", "Co": "Co6+3", "Ni": "Ni4+2", "Cu": "Cu3+1", "Zn": "Zn3+2", "Ga": "Ga3+3", "Ge": "Ge3", "As": "As3+3", "Se": "Se3+2", "Br": "Br", "Kr": "Kr4+4", "Rb": "Rb", "Sr": "Sr6+2", "Y": "Y_3+3", "Zr": "Zr3+4", "Nb": "Nb3+5", "Mo": "Mo6+6", "Tc": "Tc6+5", "Ru": "Ru6+2", "Rh": "Rh6+3", "Pd": "Pd4+2", "Ag": "Ag1+1", "Cd": "Cd3+2", "In": "In3+3", "Sn": "Sn3", "Sb": "Sb3+3", "Te": "Te3+2", "I": "I_", "Xe": "Xe4+4", "Cs": "Cs", "Ba": "Ba6+2", "La": "La3+3", "Ce": "Ce6+3", "Pr": "Pr6+3", "Nd": "Nd6+3", "Pm": "Pm6+3", "Sm": "Sm6+3", "Eu": "Eu6+3", "Gd": "Gd6+3", "Tb": "Tb6+3", "Dy": "Dy6+3", "Ho": "Ho6+3", "Er": "Er6+3", "Tm": "Tm6+3", "Yb": "Yb6+3", "Lu": "Lu6+3", "Hf": "Hf3+4", "Ta": "Ta3+5", "W": "W_6+6", "Re": "Re6+5", "Os": "Os6+6", "Ir": "Ir6+3", "Pt": "Pt4+2", "Au": "Au4+3", "Hg": "Hg1+2", "Tl": "Tl3+3", "Pb": "Pb3", "Bi": "Bi3+3", "Po": "Po3+2", "At": "At", "Rn": "Rn4+4", "Fr": "Fr", "Ra": "Ra6+2", "Ac": "Ac6+3", "Th": "Th6+4", "Pa": "Pa6+4", "U": "U_6+4", "Np": "Np6+4", "Pu": "Pu6+4", "Am": "Am6+4", "Cm": "Cm6+3", "Bk": "Bk6+3", "Cf": "Cf6+3", "Es": "Es6+3"}


atoms = read('input.pdb')

def write_itp(atoms):
    o = open('sub.itp', 'w')
    o.write('''
[ moleculetype ]
;name nrexcl
sub      3

[ atoms ]
;  nr     type resnr resnm   atom cgrp     charge      mass
''')

    for n, i in enumerate(atoms):
        q = 0.0
        atp = uff_atp[i.symbol]
        o.write('%4d%9s%5d%6s%7s%5d%12.5f%12.5f\n'%
                (n+1, atp, 1, 'sub', i.symbol, 1, q, i.mass))
    o.write('\n[ position_restraints ]\n')
    o.write('; i funct fcx fcy fcz\n')
    for n, i in enumerate(atoms):
        o.write('%d 1 100000 100000 10000\n'%(n+1))
    o.close()

write_itp(atoms)

