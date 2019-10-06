"""
@note:
"""

import os, sys

cmd = 'echo Box-X | gmx_mpi energy'
f = os.popen(cmd)
for i in f:
    if 'Box-X' in i:
        tokens = i.strip().split()
f.close()

box_x = float(tokens[1])
os.system("gmx_mpi editconf -f conf.gro -o conf.gro -box %.5f"%box_x)
