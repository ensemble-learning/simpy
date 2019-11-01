import MDAnalysis as mda

def get_elements():
    atoms_elements = []
    f = open('out.xyz', 'r')
    lines = f.readlines()
    f.close()
    for i in range(2,len(lines)):
        tokens = lines[i].strip().split()
        if len(tokens) > 3:
            ele = tokens[0]
            atoms_elements.append(ele)
    return atoms_elements
    
atoms_elements = get_elements()

tpr = 'topol.tpr'
gro = 'conf.gro'
trr = 'traj.trr'

u = mda.Universe(gro, trr)
with mda.Writer('traj.xyz', n_atoms=u.atoms.n_atoms, atoms=atoms_elements) as xyz:
    for ts in u.trajectory:
        xyz.write(ts)
