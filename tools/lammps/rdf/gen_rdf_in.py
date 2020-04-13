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

units        real
read_data       lammps.data
atom_modify     sort 0 0

pair_style      lj/cut 10
pair_coeff      * * 1 1

'''%n_bins)

for i in range(len(atomtypes)):
    o.write('group %d type %d # \n'%(i+1,i+1))
o.write('\n')

for i in range(len(atomtypes)):
    for j in range(i, len(atomtypes)):
        o.write('compute r%d%d all rdf ${n_bins} %d %d\n'%(i+1, j+1, i+1, j+1))
        o.write('fix f%d%d all ave/time 1 %d %d c_r%d%d[*] file r%s-%s.dat mode vector\n'
                %(i+1, j+1, n_frames-1, n_frames-1, i+1, j+1, atomtypes[i], atomtypes[j]))

o.write('\n')
o.write('rerun           movie.xyz dump x y z box no format xyz\n')
o.close()

