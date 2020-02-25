def read_bond_table():
    atoms = []
    f = open('test.tatb', 'r')
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split()
            atom_id = int(tokens[0])
            atom_type = int(tokens[1])
            n_bond = int(tokens[2])
            bond_atoms = [int(j) for j in tokens[3:3+n_bond]]
            atoms.append([atom_id, atom_type, n_bond, bond_atoms, -1])
    f.close()
    return atoms

def get_molecule_id(atom_id, mol_id):
    if atoms[atom_id][4] < 0:
        atoms[atom_id][4] = mol_id
        molecules[mol_id].append(atom_id)
        for partner in atoms[atom_id][3]:
            partner = partner - 1
            get_molecule_id(partner, mol_id)

atoms = read_bond_table()

molecules = []
for i in range(len(atoms)):
    molecules.append([])

for i in range(len(atoms)):
    get_molecule_id(i,i)
o = open('molecule.dat', 'w')
for mol in molecules:
    if mol:
        for i in mol:
            o.write("%d "%(i+1))
        o.write('\n')
o.close()
        
#@todo
# generate bond table
# assign elements
# assing name

