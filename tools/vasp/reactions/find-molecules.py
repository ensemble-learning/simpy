"""
@todo: assign atom types
"""
import os, sys
import numpy as np
import ase
import ase.io
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.data import covalent_radii, atomic_numbers

def reax_bond_order(a, b, r, ffield):
    bo = math.exp(pbo1*math.power(r/r0, pbo2))
    return bo

def build_cutoff_table(atps):
    scale = 1.2
    cutoff_table = {}
    for i in range(len(atps)):
        for j in range(i, len(atps)):
            ri = covalent_radii[atomic_numbers[atps[i]]]
            rj = covalent_radii[atomic_numbers[atps[j]]]
            r0 = scale*(ri+rj)
            cutoff_table[(atps[i], atps[j])] = r0
    return cutoff_table

def get_lists(atoms, cutoff_table):

    # build neighbor list
    nl = []
    for i in range(len(atoms)):
        nl.append([])

    tokens_i, tokens_j = neighbor_list('ij', atoms, cutoff_table)
    for i in range(len(tokens_i)):
        nl[tokens_i[i]].append(tokens_j[i])

    return nl

def get_molecule_id(atom_id, mol_id, atom_visits, molecules):
    if atom_visits[atom_id][2] < 0:
        atom_visits[atom_id][2] = mol_id
        molecules[mol_id].append(atom_id)
        for partner in atom_visits[atom_id][1]:
            get_molecule_id(partner, mol_id, atom_visits, molecules)

def find_molecules(atoms, nl):
    molecules = []
    for i in range(len(atoms)):
        molecules.append([])

    atom_visits = []
    for i in range(len(atoms)):
        atom_visits.append([atoms[i].symbol, nl[i], -1])

    for i in range(len(atoms)):
        get_molecule_id(i,i, atom_visits, molecules)

    molecules = [i for i in molecules if i]
    return molecules

def output_molecules(molecules, atoms):
    o = open('mol.dat', 'w')
    for i in molecules:
        for j in i:
            o.write('%5d '%j)
        o.write('\n')
    o.close()

if len(sys.argv) > 1:
    fname = sys.argv[1]
    atoms = ase.io.read(fname)
    atps = list(set(atoms.get_chemical_symbols()))
    cutoff_table = build_cutoff_table(atps)
    #cutoff_table[('O', 'Li')] = 0.0
    #cutoff_table[('Li', 'O')] = 0.0
    nl = get_lists(atoms, cutoff_table)
    molecules = find_molecules(atoms, nl)
    output_molecules(molecules, atoms)

