""" Scale the parameters in params file (ReaxFF training file)
"""

import sys
import shutil

def usage():
    print """scale the paramters in params file.
    python pscale.py [scale] [nloops]
    scale: (float) the scale parameter
    nloops: (int) number of loops
    """
def pscale(scale, nloops=1):
    lines = []
    f = open("params", "r")
    for i in f:
        if "!" in i:
            tokens = i.strip().split("!")[0].split()
            cmt = "! " + i.strip().split("!")[1] + "\n"
        else:
            tokens = i.strip().split()
            cmt = "!\n"

        if len(tokens) == 6:
            t1 = int(tokens[0])
            t2 = int(tokens[1])
            t3 = int(tokens[2])
            val = 0.5 * (float(tokens[4]) + float(tokens[5]))
            if val == 0:
                pass
            else:
                if val >0:
                    upper = val * (1+scale)
                    lower = val * (1-scale)
                    step = val * scale / 10
                else:
                    upper = val * (1-scale)
                    lower = val * (1+scale)
                    step = -val * scale / 10
                line = "%-4d%6d%6d%10.4f%10.4f%10.4f"%(t1, t2, t3, step, lower, upper)
                line += cmt
                lines.append(line)
    f.close()

    shutil.copy("params", "params.bak")
    o = open("params", "w")
    for i in range(nloops):
        for j in lines:
            o.write(j)
    o.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    elif len(sys.argv) == 2:
        print "warning: using default nloops (1)"
        scale = float(sys.argv[1])
        pscale(scale)
    elif len(sys.argv) == 3:
        scale = float(sys.argv[1])
        nloops = int(sys.argv[2])
        pscale(scale, nloops)
