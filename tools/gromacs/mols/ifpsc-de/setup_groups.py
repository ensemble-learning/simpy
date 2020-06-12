o = open('groups', 'w')
n_atoms = 1400
n_mol = 50
n_mol_atoms = int(n_atoms/n_mol)
o.write('%d %d\n'%(n_atoms, n_mol))
for i in range(n_mol):
    o.write('%d 1\n'%n_mol_atoms)
    for j in range(n_mol_atoms):
        o.write('%d '%(i*n_mol_atoms+1+j))
    o.write('\n')
o.close()
