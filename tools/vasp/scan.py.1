"""
Generate vasp isotropic scan
"""
import os
import shutil
import numpy as np
import sys

def gen_scan(start=0.8, end=1.2):
    lines = []
    f = open("POSCAR", "r")
    for i in f:
        lines.append(i) 
    f.close()
    
    start = np.power(start, 1/3.0)
    end = np.power(end, 1/3.0)
    n = 11
    
    x = np.linspace(start, end, n)
    
    for i in range(n):
        folder = "scan_%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        shutil.copy("../INCAR", ".")
        shutil.copy("../POTCAR", ".")
        shutil.copy("../KPOINTS", ".")
        shutil.copy("../pbs", ".")
        o = open("POSCAR", "w")
        lines[1] = "    %.4f\n"%x[i]
        for j in lines:
            o.write(j)
        o.close()
        os.chdir("..")

def main():
    if len(sys.argv) > 2:
        start = float(sys.argv[1])
        end = float(sys.argv[2])
        gen_scan(start, end)
    else:
        gen_scan()

if __name__ == "__main__":
    main()
