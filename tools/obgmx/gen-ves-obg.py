import os, sys
from ase.io import read, write

mid = int(sys.argv[1])
mol = sys.argv[2]
mol_name = mol.split('.')[0]

atps = []
f = open('obgmx.top', 'r')
for i in f:
    if i.strip().startswith('[ atomtypes ]'):
        break
for i in f:
    if len(i.strip()) == 0:
        break
    else:
        if not i.strip().startswith(';'):
            print(i.strip())
            atps.append(i)
f.close()

f = open('obgmx.itp', 'r')
lines = f.readlines()
f.close()

folder = '%03d-%s-opt'%(mid, mol_name)
os.makedirs(folder, exist_ok=True)

o = open('%s/%s_qforce.itp'%(folder, folder), 'w')
o.write('[ atomtypes ]\n')
for i in atps:
    o.write(i)
o.write('[ nonbond_params ]\n\n')
for i in lines:
    o.write(i)
o.close()

atoms = read(mol)
write('%s/input.pdb'%folder, atoms)
