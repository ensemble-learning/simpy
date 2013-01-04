#!/usr/bin/env python
import sys
ELEMENT = { 
        "O" : "O",
        "O1" : "O ",
        "O2" : "O ",
        "O3" : "O ",
        "Si" : "Si",
        "Si0" : "Si ",
        "Si1" : "Si ",
        "N" : "N",
        "N1" : "N ",
        "C" : "C",
        "H" : "H",
        "Na" : "Na",
        "Na0" : "Na",
        "OW" : "O",
        "HW1" : "H",
        "HW2" : "H",
        "Ar" : "Ar",
        }

def readPdb(pdbfile):
    f = open(pdbfile, 'r')
    pdbdata = []
    for i in f:
        pdbdata.append(i[:-1])
    return pdbdata

def sPdb(pdbdata):
    pdbdata_s = []
    for i in pdbdata:
        if i.split()[0] == "ATOM":
            ele = ELEMENT[i.split()[2]]
            line = i + "%12s"%ele
        else:
            line = i
        pdbdata_s.append(line)
    return pdbdata_s

def outPdb(pdbdata):
    o = open("sout.pdb", 'w')
    for j in pdbdata:
        o.write("%s\n"%j)
    o.close()

if __name__ == "__main__":
    pdbfile = sys.argv[1]
    data = readPdb(pdbfile)
    sdata = sPdb(data)
    outPdb(sdata)

