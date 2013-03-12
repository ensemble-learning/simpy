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

def main(fname, min, max, n):
    lines = []
    f = open(fname, "r")
    for i in f:
        if i.startswith("B"):
            tokens = i.strip().split()
            tokens[3] = "%bond%"
            line = " ".join(tokens)
        else:
            line = i
        lines.append(line)
    f.close()
    print lines

    bonds = np.linspace(min, max, n)

    for i in range(len(bonds)):
        o = open("scan_%02d.gjf"%i, "w")
        bond = bonds[i]
        for j in lines:
            j = j.replace("%bond%", "%.6f"%bond)
            o.write(j)
        o.close()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        usage()
    else:
        fname = sys.argv[1]
        min = float(sys.argv[2])
        max = float(sys.argv[3])
        n = int(sys.argv[4])
        main(fname, min, max, n)
