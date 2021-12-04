import os
if os.path.exists('mixture.pdb'):
    os.remove('mixture.pdb')
    
f = open('inp', 'r')
lines = f.readlines()
f.close()
mols = lines[0].strip().split()
nmols = lines[1].strip().split()
box = lines[2].strip()

o = open('pack.inp', 'w')
o.write('tolerance 2.0\n')
o.write('filetype pdb\n')
o.write('output mixture.pdb\n')
o.write('\n')
for n, i in enumerate(mols):
    o.write('structure %s\n'%i)
    o.write('    number %s\n'%nmols[n])
    o.write('    inside box 1.0 1.0 1.0 %s\n'%box)
    o.write('end structure\n\n')
o.close()

os.system('packmol < pack.inp')
lx, ly, lz = [0.1*float(i) for i in box.split()]
os.system('gmx editconf -f mixture.pdb -o input.pdb -box %.2f %.2f %.2f'%(lx, ly, lz))
os.system('gmx solvate -cp input.pdb -cs spc216')
