"""Generate trainset file
"""

import os
import sys
import numpy as np

def usage():
    print """python addtran.py [reference name] [reference energy] [energy unit]
    @note: Now we only support kcal and ev"""

if len(sys.argv) < 3:
    usage()
else:
    rname = sys.argv[1]
    rener = float(sys.argv[2])
    if len(sys.argv) == 3:
        print "Using default energy unit (kCal/mol)"
        unit = 1
    else:
        if sys.argv[3] == "ev":
            unit = 23
        else:
            print "No defined units."
            unit = 1

    ener = np.loadtxt("results")
    conf = []
    
    f = open("geo", "r")
    for i in f:
        if i.strip().startswith("DESCRP"):
            tokens = i.strip().split()
            conf.append(tokens[1])
    f.close()
    
    o = open("add.trainset", "w")
    o.write("# %s\n"%rname)
    if len(conf) == len(ener):
        for i in range(len(conf)):
            de = (ener[i] - rener) * unit 
            o.write(" 0.1 + %s/1 - %s/1  %9.2f\n"%(conf[i], rname, de))
    else:
        print "Error!"
    o.close()
    
