from ase.io import read, write

atoms = read('POSCAR')
selected = [119, 168, 49, 50, 51, 52, 53, 54]
remove = []
for i in range(len(atoms)):
    if not i+1 in selected:
        remove.append(i)

del(atoms[remove])
write('out.vasp', atoms, vasp5=True, direct=True)
