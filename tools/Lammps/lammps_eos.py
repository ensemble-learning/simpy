import shutil
import numpy as np
import os

scale = np.linspace(0.6, 1.4, 11)
scale = np.power(scale, 1/3.0)

line = "change_box all x scale %sc% y scale  %sc% z scale %sc% xy scale %sc% xz scale %sc% yz scale %sc% remap\n"

lines = []
f = open("lammps_input", "r")

counter = 0
for i in f:
    if i.strip().startswith("#"):
        pass
    elif "relax" in i:
        pass
    else:
        lines.append(i)
        if "read_restart" in i or "read_data" in i:
            lines.append(line)
            n = counter + 1
        counter += 1

for i in range(len(scale)):
    folder = "eos_%02d"%i
    if not os.path.exists(folder):
        os.mkdir(folder) 
    os.chdir(folder)
    sc = scale[i]
    print sc    
    lines[n] = line.replace("%sc%", "%.3f"%sc)
    o = open("lammps_input", "w")
    for j in lines:
        o.write(j)
    o.close()
    shutil.copy("../min.rst", ".")
    shutil.copy("../ffield", ".")
    os.chdir("..")
