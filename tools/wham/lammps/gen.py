import os, shutil, time
import numpy as np

plumed = """
UNITS LENGTH=A ENERGY=kcal/mol
dd: DISTANCE ATOMS=1,2
RESTRAINT ARG=dd KAPPA=%force% AT=%distance%
PRINT ARG=dd FILE=colvar STRIDE=100
"""

force = 25
data = np.linspace(2.5,6, 18)

for i in range(len(data)):
    folder = 'r%02d'%i
    if not os.path.exists(folder):
        os.mkdir(folder)
    if 0:
        shutil.copy('in.lmp', folder)
        if i == 0:
            shutil.copy('data.lmp.0', './%s/data.lmp'%folder)
        else:
            shutil.copy('data.lmp', folder)
    os.chdir(folder)
    o = open('plumed.dat', 'w')
    control = plumed.replace('%distance%', '%.4f'%data[i])
    control = control.replace('%force%', '%.4f'%force)
    o.write(control)
    o.close()
    if 0:
        os.system('/central/home/tcheng/soft/lammps/lammps-3Mar20/src/lmp_mpi -in in.lmp')
        time.sleep(1)
        shutil.copy('data.eq.lmp', '../data.lmp')
    os.chdir('../')

o = open('metadata', 'w')
for i in range(len(data)):
    folder = 'r%02d'%i
    o.write('./%s/colvar '%folder)
    o.write('%.4f %.2f \n'%(data[i], 2*force))
o.close()

o = open('ana.sh', 'w')
o.write('wham %.2f %.2f 100 0.01 300 0 metadata out1.dat 1000 1\n'%(data[0], data[-1]))
o.write('wham %.2f %.2f 50 0.01 300 0 metadata out2.dat 1000 1\n'%(data[0], data[-1]))
o.write('wham %.2f %.2f 25 0.01 300 0 metadata out3.dat 1000 1\n'%(data[0], data[-1]))
o.close()

