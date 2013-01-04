#!/usr/bin/env python

import os
import math
import sys

def crystal(folder):
    a = []
    f1 = None
    f2 = None
    f3 = None
    f1, f2, f3 = os.popen3(" echo Box-XX Box-YX Box-YY Box-ZX Box-ZY Box-ZZ| g_energy_mpi -f %s"%folder)
    for i in f2 :
        if i[:3] == "Box" :
            a.append(float(i.split()[1]))
    f2.close()
    print folder

    box_xx = abs( a[0] )
    box_yy = math.sqrt(a[1]*a[1] + a[2]*a[2])
    box_zz = math.sqrt(a[3]*a[3] + a[4]*a[4] + a[5]*a[5])
    alpha = math.acos((a[1]*a[3] + a[2]*a[4]) / box_xx / box_zz) / math.pi * 180
    beta = math.acos(a[0]*a[3] / box_xx / box_zz) / math.pi * 180
    gama = math.acos(a[0]*a[1] / box_xx / box_yy) / math.pi * 180
    return [box_xx,box_yy,box_zz,alpha,beta,gama]

def parse(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if i.endswith(".edr"):
            cellParameter[os.path.basename(folder)] = crystal(fullname)
        elif os.path.isdir(fullname):
            parse(fullname)

def output(dict):
    from cell_parameter import EXP
    keys = dict.keys()
    keys.sort()
    o = open("cell.log", 'w')
    for i in keys:
        for j in range(len(dict[i])):
            o.write("%20s"%i)
            o.write("%9.3f%9.3f"%(dict[i][j], EXP[i][j]))
            o.write("%9.3f\n"%(dict[i][j]- EXP[i][j]))
    o.close()
if __name__ == "__main__":
    global cellParameter
    cellParameter = {}
    parse(".")
    output(cellParameter)
