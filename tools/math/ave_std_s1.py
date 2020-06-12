import numpy as np
import sys

if len(sys.argv) > 1:
    data = np.loadtxt(sys.argv[1])
    print(np.average(data))
    print(np.std(data))
