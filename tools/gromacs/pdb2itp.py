"""
@todo: assign atom types
"""
import os, sys
import numpy as np
import ase
import ase.io
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.geometry import distance, get_angles
from ase.data import covalent_radii, atomic_numbers

cutoff_table = {('H', 'H'):1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85, ('H', 'N'):1.3, ('N', 'N'):1.85, ('C', 'N'):1.85}
#atom_types_table = {'C': 'opls_135', 'H': 'opls_140', 'N': 'opls_237', 'Cs': 'Cs', 'Sn': 'Sn', 'I':'I'}
atom_types_table = {'C': 'C', 'H': 'H', 'N': 'N', 'Cs': 'Cs', 'Sn': 'Sn', 'I':'I', 'Au':'Au', 'S':'S', 'O': 'O'}
ERROR_ATP = 'Warning: User-defined atom types provided but the numbers do not match!\n'

class Params():
    def __init__(self,):
        self.fname = ''
        self.cutoff_file = ''
        self.cutoff_table = {}
        self.nl = []
        self.molecules = []
        self.tolerance = 1.1
        self.ignore_element = []

def update_cutoff_table(cutoff_table):
    if os.path.exists(p.cutoff_file):
        f = open(p.cutoff_file, 'r')
        for i in f:
            tokens = i.strip().split()
            if len(tokens) == 3:
                p.cutoff_table[(tokens[0], tokens[1])] = float(tokens[2])
                p.cutoff_table[(tokens[1], tokens[0])] = float(tokens[2])
        f.close()

def build_cutoff_table(atps, p):
    for i in range(len(atps)):
        for j in range(i, len(atps)):
            if atps[i] in p.ignore_element or atps[j] in p.ignore_element:
                r0 = 0.1
            else:
                ri = covalent_radii[atomic_numbers[atps[i]]]
                rj = covalent_radii[atomic_numbers[atps[j]]]
                r0 = p.tolerance*(ri+rj)
            p.cutoff_table[(atps[i], atps[j])] = r0

def get_lists(atoms, p):

    # build neighbor list
    nl = []
    for i in range(len(atoms)):
        nl.append([])

    tokens_i, tokens_j = neighbor_list('ij', atoms, p.cutoff_table)
    for i in range(len(tokens_i)):
        nl[tokens_i[i]].append(tokens_j[i])

    o = open('bonds.dat', 'w')
    # build bond list
    bond_list = []
    for i in range(len(nl)):
        if len(nl[i]) > 0:
            ai = i
            for j in range(len(nl[i])):
                aj = nl[i][j]
                if ai < aj:
                    bond_length = atoms.get_distance(ai, aj, mic=True)
                    bond_length = bond_length/10.0
                    bond_list.append([ai, aj, bond_length])
                    o.write('%d %d %.4f\n'%(ai, aj, bond_length))
    o.close()
    n_bond = len(bond_list)
    print(n_bond, "bond terms")

    o = open('angles.dat', 'w')
    angle_list = []
    # build angle list
    for i in range(len(nl)):
        if len(nl[i]) > 1:
            aj = i
            for j in range(len(nl[i])):
                ai = nl[i][j]
                for k in nl[i][j+1:]:
                    ak = k
                    angle = atoms.get_angle(ai, aj, ak, mic=True)
                    angle_list.append([ai, aj, ak, angle])
                    o.write('%d %d %d %.4f\n'%(ai, aj, ak, angle))
    o.close()
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
                        angle = atoms.get_dihedral(di, dj, dk, dl, mic=True)
                        dihedral_list.append([di, dj, dk, dl, angle])
    n_dihedral = len(dihedral_list)
    print(n_dihedral, "diheral terms")
    return nl, bond_list, angle_list, dihedral_list

def to_itp(fname, atoms, nl, bl, al, dl, charges, atomtypes):
    o = open(fname, "w")
    molname = fname.split(".")[0]
    chemical_symbols = atoms.get_chemical_symbols()
    masses = atoms.get_masses()
    o.write("\n[ moleculetype ]\n")
    o.write("; Name    nrexcl\n")
    o.write(" %s    3\n"%molname)
    o.write("\n[ atoms ]\n")
    o.write("; nr  type  resnr  residue  atom  cgnr  charge  mass\n")
    for i in range(len(atoms)):
        o.write("%6d    %s    1    %s"%((i+1), atomtypes[i], molname))
        o.write("%4s    1    %10.6f %10.4f\n"%(chemical_symbols[i], charges[i], masses[i]))

    o.write("\n[ bonds ]\n")
    o.write("; ai  aj  funct  c0  c1  c2  c3\n")
    for i in bl:
        o.write("%6d%6d    1"%(i[0]+1, i[1]+1))
        o.write(' %12.6f %12.2f\n'%(i[2], 270000)) 

    o.write("\n[ pairs ]\n")
    o.write("; ai  aj  ak  al  funct  phi  cp mult\n")
    for i in dl:
        o.write("%6d%6d\n"%(i[0]+1, i[3]+1))
    
    o.write("\n[ angles ]\n")
    o.write("; ai  aj  ak  funct  c0  c1  c2  c3\n")
    for i in al:
        o.write("%6d%6d%6d    2"%(i[0]+1, i[1]+1, i[2]+1))
        o.write(' %8.2f %12.2f\n'%(i[3], 1200)) 

    o.write("\n[ dihedrals ]\n")
    o.write("; ai  aj  ak  al  funct  phi  cp  mult\n")
    for i in dl:
        o.write("%6d%6d%6d%6d    2"%(i[0]+1, i[1]+1, i[2]+1, i[3]+1))
        o.write(' %8.2f %12.2f\n'%(i[4], 0.1)) 
    o.write("\n[ system ]\n")
    o.write(";name\n%s\n\n"%molname) 
    o.write("[ molecules ]\n\n")
    o.write("; Compound    #mols\n") 
    o.write("%s    0\n"%molname)
    o.close()

def assign_charges(atoms):
    charges = [0.0]*len(atoms)
    if os.path.exists("q.dat"):
        charges = np.loadtxt("q.dat")
    #print("%.4f"%np.sum(charges), "total charges")
    return charges   

def assign_atom_types_opls(atoms, nl):
    chemical_symbols = atoms.get_chemical_symbols()
    atomtypes, atomtypes_usr = [], []
    for i in range(len(atoms)):
        atomtypes.append(atom_types_table[chemical_symbols[i]])
    if os.path.exists('atp.dat'):
        f = open('atp.dat', 'r')
        for i in f:
            tokens = i.strip()
            if len(tokens) > 0:
                atomtypes_usr.append(tokens)
        f.close()
        if len(atomtypes_usr) == len(atomtypes):
            for j in range(len(atomtypes)):
                atomtypes[j] = atomtypes_usr[j]
        else:
            sys.stderr.write(ERROR_ATP)
        
    return atomtypes

if len(sys.argv) > 1:
    fname = sys.argv[1]
    itpfile = fname.split(".")[0] + ".itp"

    p = Params()

    atoms = ase.io.read(fname)
    atps = list(set(atoms.get_chemical_symbols()))
    build_cutoff_table(atps, p)
    update_cutoff_table(p)

    nl, bl, al, dl = get_lists(atoms, p)
    charges = assign_charges(atoms)
    atomtypes = assign_atom_types_opls(atoms, nl)

    to_itp(itpfile, atoms, nl, bl, al, dl, charges, atomtypes)

