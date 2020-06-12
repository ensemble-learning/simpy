import numpy as np
import sys

fname = sys.argv[1]

data = np.loadtxt(fname)
data = data.transpose()

for i in data:
    print(np.average(i))
