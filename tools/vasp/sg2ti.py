#!/usr/bin/env python

import os, shutil, time
def zero_increm():
    stamp = int(time.time()*100)
    stamp = str(stamp)

    f = open("INCAR", "r")
    lines = f.readlines()
    f.close()

    for i in range(len(lines)):
        if lines[i].strip().startswith("INCREM"):
            lines[i] = "INCREM    =    0.0\n"
    shutil.copy("INCAR", "INCAR.%s"%stamp[-8:])
    o = open("INCAR", "w")
    for i in lines:
        o.write(i)
    o.close()

os.system("xdat2pos")
os.chdir("./poscars")
os.system("select_frame 1 2000 11")
os.chdir("wf")
os.system("ls POSCAR_000* > flist")
os.system("bsr")
os.system("rm POSCAR_* flist")
os.system("cmdall 'cp POSCAR_* POSCAR'")
shutil.copy("../../KPOINTS", ".")
shutil.copy("../../INCAR", ".")
shutil.copy("../../POTCAR", ".")
shutil.copy("../../ICONST", ".")
shutil.copy("../../pbs", ".")
zero_increm()
os.system("copyall KPOINTS INCAR POTCAR ICONST pbs")
os.chdir("..")
os.system("mv wf ../../ti")


