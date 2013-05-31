"""
get energy terms from fort.73
"""

import os
import numpy as np
import matplotlib.pyplot as plt

assert os.path.exists("fort.73")

lines = []
f = open("fort.73", "r")
for i in f:
    lines.append(i)
f.close()

data = []
for i in range(len(lines)/3):
    tokens = lines[3*i+2].strip().split()
    data.append([float(j) for j in tokens])

data = np.array(data)
print data
    

