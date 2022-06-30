import os
'''
obabel -:"CCCCCCCCCCC(C(CCCCC)CCCCCC)C" -omol2  -h --gen3d > out.mol2
mol.graph.edges[0, 1]['order'] = 2
mol2 file format has bond order information
'''

import copy
import selfies as sf
import networkx as nx
from pysmiles import write_smiles, fill_valence, remove_explicit_hydrogens
import matplotlib.pyplot as plt

class Molecule():
    def __init__(self,):
        self.chemical_symbols = []
        self.bonds = []
        self.graph = nx.Graph()

def write_mols(mol, fname='out.pdb'):

    largest_cc = max(nx.connected_components(mol), key=len)

    mol = mol.subgraph(largest_cc).copy()
    nx.draw(mol, with_labels = True)
    #plt.show()

    mol_smi_no_h = write_smiles(mol)
    fill_valence(mol, respect_bond_order=True, respect_hcount=True, max_bond_order=4)
    #fill_valence(mol, respect_hcount=True)
    mol_smi = write_smiles(mol)

    try:
        mol_sf = sf.encoder(mol_smi_no_h)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]
        sf_smi = sf.decoder(mol_sf)  # C1=CC=CC=C1
        os.system('obabel -:"%s" -opdb  -h --gen3d > %s'%(mol_smi, fname))
    except sf.EncoderError:
        pass  # sf.encoder error!
    except sf.DecoderError:
        pass  # sf.decoder error!

    print(mol_smi)

mol = Molecule()
f = open('atoms', 'r')
for i in f:
    tokens = i.strip().split()
    if len(tokens) > 0:
        mol.chemical_symbols.append(tokens[1])
f.close()

f = open('bonds', 'r')
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 2:
        mol.bonds.append((int(tokens[0]), int(tokens[1]), 1))
    elif len(tokens) == 3:
        mol.bonds.append((int(tokens[0]), int(tokens[1]), int(tokens[2])))
f.close()

mol2 = nx.Graph()
for (i,j,k) in mol.bonds: 
    mol2.add_edge(i-1,j-1)
    mol2.edges[i-1,j-1]['order'] = k

for idx in range(len(mol.chemical_symbols)):
    mol2.nodes[idx]['element'] = mol.chemical_symbols[idx]

write_mols(mol2)
