"""generate the scan files for g03 calculation
"""

import os
import shutil
import sys
import numpy as np
import argparse

def angle_scan(fname, opt, min, max, n):
    lines = []
    f = open(fname, "r")
    for i in f:
        if i.startswith("A"):
            tokens = i.strip().split()
            tokens[4] = "%angle%"
            line = " ".join(tokens)
        else:
            line = i.strip() + '\n'
        lines.append(line)
    f.close()

    angles = []
    flag = 0
    for i in np.linspace(min, max, n):
        if opt < i and flag == 0:
            angles.append(opt)
            flag = 1
        angles.append(i)
    print angles
    
    oangle = open("angles", "w")
    for i in range(len(angles)):
        o = open("scan_%02d.gjf"%i, "w")
        angle = angles[i]
        for j in lines:
            j = j.replace("%angle%", "%.6f"%angle)
            o.write(j)
        o.write("\n")
        o.write("\n")
        o.close()
        oangle.write("%12.6f\n"%angle)
        if not os.path.exists("scan_%02d"%i):
            os.mkdir("scan_%02d"%i)
        shutil.copy("scan_%02d.gjf"%i, "scan_%02d"%i)
    oangle.close()

def bond_scan(fname, opt, min, max, n):
    lines = []
    f = open(fname, "r")
    for i in f:
        if i.startswith("B"):
            tokens = i.strip().split()
            tokens[3] = "%bond%"
            line = " ".join(tokens)
        else:
            line = i.strip() + '\n'
        lines.append(line)
    f.close()

    bonds = []
    start = opt - 0.4
    end = opt + 0.4
    if min < start:
        for i in np.linspace(min, start, n):
            bonds.append(i)
    for i in np.linspace(start + 0.1, end - 0.1, 7):
        bonds.append(i)
    if max > end:
        for i in np.linspace(end, max, n):
            bonds.append(i)

    obond = open("bonds", "w")
    for i in range(len(bonds)):
        o = open("scan_%02d.gjf"%i, "w")
        bond = bonds[i]
        for j in lines:
            j = j.replace("%bond%", "%.6f"%bond)
            o.write(j)
        o.write("\n")
        o.write("\n")
        o.close()
        obond.write("%12.6f\n"%bond)
        if not os.path.exists("scan_%02d"%i):
            os.mkdir("scan_%02d"%i)
        shutil.copy("scan_%02d.gjf"%i, "scan_%02d"%i)
    obond.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", nargs=1, help="g03 file")
    parser.add_argument("-params", nargs="+", type=float, help="Parameters for scan: opt, min, max, n")
    parser.add_argument("-type", default="bond", nargs=1, help="Scan type: bond or angle")
    args = parser.parse_args()
    
    fname = args.fname[0]
    print "Processing %s"%fname
    if args.params:
        tokens = args.params
        if len(tokens) < 3:
            print "Need at least 3 paramester. Only %d input."%len(tokens)
            exit()
        else:
            opt = tokens[0]
            min = tokens[1]
            max = tokens[2]
            print "opt = %.2f, min = %.2f, max = %.2f"%(opt, min, max),
            if len(tokens) > 3:
                n = tokens[3]
                print "n = %d"%n
            else:
                print ""

    if args.type:
        type = args.type[0]
    print "%s scan"%type
    print n
    if type == "bond": 
        bond_scan(fname,opt, min, max, n)
    elif type == "angle":
        angle_scan(fname, opt, min, max, n)

if __name__ == "__main__":
    main()
