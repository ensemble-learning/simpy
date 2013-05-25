"""generate the scan files for g03 calculation
"""

import sys
import numpy as np

def usage():
    print """scan.py inputfile min max n
    inputfile: model file
    min: the min value
    max: the max value
    n: number of data points
    """

def main(fname, opt, min, max, n):
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
    start = opt - 0.3
    end = opt + 0.3
    if min < start:
        for i in np.linspace(min, start, n):
            bonds.append(i)
    for i in np.linspace(start, end, 7):
        bonds.append(i)
    if max > end:
        for i in np.linspace(end, max, n):
            bonds.append(i)

    for i in range(len(bonds)):
        o = open("scan_%02d.gjf"%i, "w")
        bond = bonds[i]
        for j in lines:
            j = j.replace("%bond%", "%.6f"%bond)
            o.write(j)
        o.write("\n")
        o.write("\n")
        o.close()


if __name__ == "__main__":
    if len(sys.argv) < 5:
        usage()
    else:
        fname = sys.argv[1]
        opt = float(sys.argv[2])
        min = float(sys.argv[3])
        max = float(sys.argv[4])
        if len(sys.argv) == 6:
            n = int(sys.argv[4])
        else:
            n = 3
        main(fname,opt, min, max, n)
