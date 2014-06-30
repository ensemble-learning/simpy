"""
Calculate the mismatch of two surface
"""

import sys

def cal_mismatch(sa, sb, la, lb):
    lines = []
    assert sa < la, "sa %.4f should smaller than la %.4f"%(sa, la)
    assert sb < lb, "sb %.4f should smaller than lb %.4f"%(sb, lb)
    #print sa, la, sb, lb
    nmax = 10
    for i in range(1,nmax):
        for j in range(1,nmax):
            anew = la*i
            bnew = lb*j
            na = int(anew/sa)
            #print anew, na
            a = na * sa
            nb = int(bnew/sb)
            #print bnew, nb
            b = nb * sb
            dev = (1 - sa*sb*na*nb/(anew*bnew))*100
            #print na, nb, anew, bnew, sa, la, sb, lb
            line = [na, nb, dev, anew, bnew, i, j]
            lines.append(line)
    data = sorted(lines, key=lambda tmp:tmp[2])
    for i in data:
        print "%8.2f%8d%8d%8d%8d%8.2f%8.2f"%(i[2], i[0], i[1], i[5], i[6], i[3], i[4])

def main():
    if len(sys.argv) < 5:
        print __doc__
    else:
        sa = float(sys.argv[1])
        sb = float(sys.argv[2])
        la = float(sys.argv[3])
        lb = float(sys.argv[4])
        cal_mismatch(sa, sb, la, lb)

if __name__ == "__main__":
    main()
