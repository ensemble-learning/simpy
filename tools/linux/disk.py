import os
import time

FOLDERS = """ADF
ASE
BENCHMARK
bin
bio
BUILD
CCDC
CIF
CRYSTAL
DATABASE
DEBUG
DFFJOBS
DFTB
GROMACS
GROUP
GULP
JAGUAR
LAMMPS
lib
Mail
PUREMD
QM
reaxFF
SHARE
soft
TMP
TPS
TRANSFER
VAR
VASP
WORKSPACE
"""

DEBUG = 0

def get_size():
    data = []
    for i in FOLDERS.split():
        os.chdir(i)
        for j in os.popen("du -s --block-size=1M"):
            tmp = "%09d_%s"%(int(j.strip().split()[0]), i)
            data.append(tmp)
        os.chdir("..")
    data.sort(reverse=True)
  
    f = open("/net/hulk/home6/chengtao/disk.log", "w")
    for i in data:
        tokens = i.strip().split("_")
        size = int(tokens[0])
        name = tokens[1]
        if DEBUG:
            print "%20s%20d"%(name, size)
        f.write("%20s%20d\n"%(name, size))
    f.close()
    os.system("chmod g-r /net/hulk/home6/chengtao/disk.log")

if __name__ == "__main__":
    while(1):
        get_size()
        time.sleep(21600)
