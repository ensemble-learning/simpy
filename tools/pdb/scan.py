"""
Generate input scan file of G03
"""
import numpy as np
import os
import sys

mol = sys.argv[1]

a = np.linspace(0.5, 1.0, 11)
b = np.linspace(1.0, 2.0, 11)
data = []

folder = "%s_cal"%mol

if not os.path.exists(folder):
    os.mkdir(folder)
for i in a:
    data.append(i)

for i in b[1:]:
    data.append(i)

for i in range(len(data)):
    scale = i
    os.system("editconf -f %s.pdb -o scan_%02d.pdb -scale %.4f"%(mol, i, data[i]))
    os.system("babel -ipdb scan_%02d.pdb -ogjf scan_%02d.gjf"%(i, i))
    f = open("scan_%02d.gjf"%(i))
    lines = f.readlines()
    lines[0] = "# b3lyp/6-311g(d,p)\n"
    lines[4] = " 0 1\n"
    outfile = os.path.join(os.getcwd(), folder)
    outfile = os.path.join(outfile, "scan_%02d_cal.gjf"%i)
    o = open(outfile, "w")
    for j in lines:
        o.write(j)
    o.close()
    f.close()

