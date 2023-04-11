from ase.io import read, write

class Mol():
    def __init__(self,):
        self.name = 0
        self.nmol = 0
        self.natoms = 0

def gen_masses(atoms):
    elements, masses = [], []
    for i in atoms:
        if i.symbol not in elements:
            elements.append(i.symbol)
            masses.append(i.mass)

    o = open('masses', 'w')
    for n, i in enumerate(elements):
        o.write('%s %.4f\n'%(elements[n], masses[n]))
    o.close()

def read_top():
    mols = {}
    mols_list = []
    itps = []
    f = open('topol.top', 'r')
    for i in f:
        if i.startswith('#include'):
            tokens = i.strip().split()
            itp_file = tokens[1][1:-1]
            itps.append(itp_file)
    f.close()

    f = open('topol.top', 'r')
    for i in f:
        if '[ molecules ]' in i:
            break
    for i in f:
        if not i.strip().startswith(';'):
            tokens = i.strip().split()
            if len(tokens) > 1:
                mol = Mol()
                mol.name = tokens[0].lower()
                mol.nmol = int(tokens[1])
                mols[mol.name] = mol
                mols_list.append(mol.name)
    f.close()
    print(mols_list)

    for i in itps:
        mol_name = ''
        natoms = 0
        f = open(i, 'r')
        for i in f:
            if i.startswith('[ moleculetype ]'):
                break
        for i in f:
            if i.startswith('[ atoms ]'):
                break
            else:
                if not i.startswith(';'):
                    tokens = i.strip().split()
                    if len(tokens) == 2:
                        mol_name = tokens[0].lower()
        for i in f:
            if i.startswith('['):
                break
            if not i.startswith(';'):
                tokens = i.strip().split()
                if len(tokens) >= 7:
                    natoms += 1
        if mol_name in mols.keys():
            mols[mol_name].natoms = natoms
        f.close()
    return mols, mols_list

def gen_groups(atoms, mols, mols_list):
    o = open('groups', 'w')
    o.write('%d %d\n'%(len(atoms), len(mols_list)))
    nn = 1
    for n, i in enumerate(mols_list):
        na = mols[i].natoms * mols[i].nmol
        o.write('%d 1\n'%na)
        for n in range(nn, na+nn):
            o.write('%d '%n)
        o.write('\n')
        nn += na
    o.close()

def gen_supergroups(atoms, mols, mols_list):
    o = open('supergroups', 'w')
    o.write('1-%d\n'%len(mols_list))
    o.close()

atoms = read('conf.gro')
gen_masses(atoms)
mols, mols_list = read_top()
gen_groups(atoms, mols, mols_list)
gen_supergroups(atoms, mols, mols_list)

