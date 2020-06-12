import os
o = open('PRECALC', 'w')
o.write("""MINDISTANCE=24.99999999872511509929
INCLUDEGAMMA=TRUE
HEADER=Verbose
""")
o.close()

os.system('demo_kplib POSCAR PRECALC > KPOINTS')


