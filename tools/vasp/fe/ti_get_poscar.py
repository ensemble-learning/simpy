import numpy as np
import shutil
import os

data = np.loadtxt("../lammps/windows.dat")
for i in range(len(data)):
    n = int(data[i][0]) + 1
    pos = "POSCAR_%06d"%n
    folder = "run_%02d"%i
    if not os.path.exists(folder):
        os.mkdir(folder)
    shutil.copy("../poscars/%s"%pos, folder)

