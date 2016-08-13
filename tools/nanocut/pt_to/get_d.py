import numpy as np
import math
import sys

if len(sys.argv) < 2:
    pass
else:
    xyzfile = sys.argv[1]
    f = open(xyzfile, "r")
    lines = f.readlines()
    coords = lines[2:]
    data = []
    for i in coords:
        data.append([float(j) for j in i.strip().split()[1:]])
        
    data = np.array(data)
    data = data.transpose()
    com = [np.average(data[0]), np.average(data[1]),
            np.average(data[2])]
    rmax = 0
    for i in data:
        r2 = i[0]*i[0] + i[1]*i[1] + i[2]*i[2]
        r = math.sqrt(r2)
        if r > rmax:
            rmax = r
    dmax = 2 * rmax
    print dmax
