!/usr/bin/env python

import os, shutil

n = 0
for i in os.listdir("."):
    if os.path.isdir(i):
        os.chdir(i)
        poscars = []
        for j in os.listdir("."):
            if j.startswith("POSCAR"):
                poscars.append(j)
        poscars.sort()
        for k in poscars:
            shutil.copy(k, "../POSCAR_%06d"%n)
            n += 1
        os.chdir("..")

        
