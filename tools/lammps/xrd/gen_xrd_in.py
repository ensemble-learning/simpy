import os
n_bins = 100
# determine atom types
atomtypes = []
f = open('lammps.data', 'r')
for i in f:
    if i.strip().startswith('Masses'):
        break

for i in f:
    if i.strip().startswith('Atoms'):
        break
    else:
        tokens = i.strip().split()
        if len(tokens) > 0:
            atomtypes.append(tokens[-1])

# determine the number of frames
f = open('movie.xyz', 'r')
lines = f.readlines()
f.close()
n_atoms = int(lines[0])
n_frames = len(lines)/(n_atoms+2)

o = open('lammps.in', 'w')
o.write('''variable n_bins equal %d

atom_style      charge
boundary        p p p

units       metal
timestep    0.001

read_data       lammps.data
atom_modify     sort 0 0

pair_style      none
atom_modify     sort 0 0

'''%n_bins)

o.write('compute XRD all xrd 1.541838')
for i in range(len(atomtypes)):
    o.write(' %s'%atomtypes[i])
o.write(' 2Theta 5 50 c 1 1 1 LP 1 echo\n')
o.write('fix             1 all ave/histo/weight 1 %d %d 5 50 500 c_XRD[1] c_XRD[2] &\n'%(n_frames-1, n_frames-1))
o.write('                mode vector file xrd.dat\n')

o.write('rerun           movie.xyz dump x y z box no format xyz\n')
o.close()

