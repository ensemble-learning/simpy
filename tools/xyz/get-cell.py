import os, sys
import numpy as np

f = open(sys.argv[1], 'r')
lines = f.readlines()
f.close()

data = []
for i in lines[2:]:
    tokens = i.strip().split()
    if len(tokens) > 3:
        x = float(tokens[1])
        y = float(tokens[2])
        z = float(tokens[3])
        data.append([x, y, z])

data = np.array(data)
data = data.transpose()

dx = np.max(data[0]) -np.min(data[0])
dy = np.max(data[1]) -np.min(data[1])
dz = np.max(data[2]) -np.min(data[2])

"""
dx = dx/10
dy = dy/10
dz = dz/10
"""

dx = dx + 40
dy = dy + 40
dz = 30
print("python ~/soft/simpy/lib/e_2_xyz.py -f *.xyz -pbc %.2f %.2f %.2f 90 90 90 -c"%(dx, dy, dz))
os.system("python ~/soft/simpy/lib/e_2_xyz.py -f *.xyz -pbc %.2f %.2f %.2f 90 90 90 -c"%(dx, dy, dz))
