""" Add keywords into geo file
"""

import shutil
import sys

def usage():
    print """python addkeys.py keyword atom1 atom2 start step
    usage: add keys in the geo file (Now only support bond restraint.
    keyword = "bond" or "angle"
    """

def addBondRes(a1, a2, start, step):

    lines = []
    com1 = "#              At1 At2   R12    Force1  Force2  dR12/dIter(MD) Start (MD) End (M\n"
    bond = 0
    counter = 0

    f = open("geo", "r")

    for i in f:
        lines.append(i)
        if i.strip().startswith("REMARK"):
            bond = start + counter * step
            key1 ="BOND RESTRAINT    %d   %d  %.4f 7500.00  1.0000  0.0000000       0       0\n"%(a1, a2, bond)
            lines.append(key1)
            counter += 1
        elif i.strip().startswith("BOND RESTRAINT"):
            print "BOND RESTRAINTs have been added"
            return 0
    f.close()

    shutil.copy("geo", "geo.addkey.bak")

    o = open("geo", "w")
    for i in lines:
        o.write(i)
    o.close()


def addAngRes(a1, a2, a3, start, step):

    lines = []
    com1 = "#              At1 At2 At3   R12    Force1  Force2  dR12/dIter(MD) Start (MD) End (M\n"
    ang = 0
    counter = 0

    f = open("geo", "r")

    for i in f:
        lines.append(i)
        if i.strip().startswith("REMARK"):
            ang = start + counter * step
            key1 = "ANGLE RESTRAINT    %d   %d   %d  %.2f  500.00 10.0000 0.000000\n"%(a1, a2, a3, ang)
            lines.append(key1)
            counter += 1
        elif i.strip().startswith("ANGLE RESTRAINT"):
            print "ANGLE RESTRAINTs have been added"
            return 0
    f.close()

    shutil.copy("geo", "geo.addAng.bak")

    o = open("geo", "w")
    for i in lines:
        o.write(i)
    o.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    else:
        key = sys.argv[1]
        if key == "bond":
            if len(sys.argv) < 6:
                print "error bond input"
            else:
                a1 = int(sys.argv[2])
                a2 = int(sys.argv[3])
                start = float(sys.argv[4])
                step = float(sys.argv[5])
                addBondRes(a1, a2, start, step)
        elif key == "angle":
            if len(sys.argv) < 6:
                print "error bond input"
            else:
                a1 = int(sys.argv[2])
                a2 = int(sys.argv[3])
                a3 = int(sys.argv[4])
                start = float(sys.argv[5])
                step = float(sys.argv[6])
                addAngRes(a1, a2, a3, start, step)

