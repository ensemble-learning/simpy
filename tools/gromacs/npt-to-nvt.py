import os
os.system('echo Box-X | gmx energy | grep Box-X > box-x.dat')

if os.path.exists('box-x.dat'):
    f = open('box-x.dat', 'r')
    lines = f.readlines()
    f.close()

tokens = lines[0].strip().split()
os.system('gmx editconf -f confout.gro -o conf.gro -box %s'%tokens[1])

