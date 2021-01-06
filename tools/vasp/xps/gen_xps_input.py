import os, shutil

o = open('KPOINTS', 'w')
o.write('''K-Points
 0
Monkhorst Pack
1  1  1
 0  0  0
''')
o.close()

folder = 'poscars'
flist = []
for i in os.listdir(folder):
    if i.strip().startswith('POSCAR_'):
        flist.append(i)

for i in flist:
    os.makedirs(i, exist_ok=True)
    shutil.copy(os.path.join(folder, i), os.path.join(i, 'POSCAR'))
    os.chdir(i)
    if not os.path.exists('KPOINTS'):
        os.system('ln -s ../KPOINTS .')
    if not os.path.exists('POTCAR'):
        os.system('ln -s ../POTCAR .')
    os.chdir('..')

    
