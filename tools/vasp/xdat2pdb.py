from ase.io import read, write

atoms = read('../POSCAR')
write('input.pdb', atoms)

traj = read('../XDATCAR', index='')
write('output.pdb', traj)

