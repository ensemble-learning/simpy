from ase.io import read, write
import os, random, copy
import numpy as np

insert = read('bf3.pdb')
com = insert.get_center_of_mass()
insert.translate(-com)

atoms = read('pt-nw-new.pdb')
cell = atoms.get_cell()
pbc = []
ngx, ngy = 8, 8
box = np.zeros([ngx, ngy])

for i in range(3):
    pbc.append(cell[i][i])
[lx, ly, lz] = pbc
dx, dy = lx/ngx, ly/ngy

for n, i in enumerate(atoms):
    [x, y, z] = i.position
    nx = int(x/dx)
    ny = int(y/dy)
    box[nx][ny] = 1

n_box_occ = np.sum(box)
nmol = ngx*ngy-n_box_occ

for i in range(ngx):
    for j in range(ngy):
        if box[i][j] == 0:
            x = i*dx+0.5*dx
            y = j*dy+0.5*dy
            z = random.random()*lz
            new = copy.deepcopy(insert)
            new.translate([x, y, z])
            atoms.extend(new)
write('p1.pdb', atoms)
os.system('gmx insert-molecules -f p1.pdb -ci clo4.pdb -nmol %d -o p2.pdb'%nmol)
os.system('gmx solvate -cp p2.pdb -cs spc216 -o p3.pdb')

