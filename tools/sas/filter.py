"""
"""
import numpy as np
import sys

cut_off = 0.04
if len(sys.argv) > 1:
    cut_off = float(sys.argv[1])
data = np.loadtxt("frac_atom.dat")
sur = np.where(data > cut_off, 1, 0)
n_sur = np.sum(sur)
n_total = len(data)
print n_sur, n_total
print n_sur*100.0/n_total
np.savetxt("sur_sas.dat", sur, fmt="%.6f")
