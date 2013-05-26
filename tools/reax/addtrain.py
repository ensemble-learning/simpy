"""Generate trainset file
"""

import os
import sys
import numpy as np

def usage():
    print
    print """python addtran.py [reference name] [reference energy] [energy unit]
    @note: Now we support kcal, hartree and ev"""

def preprocess():
    print
    print "-----------------------hint----------------------"
    flag = 1
    if os.path.exists("results"):
        data = np.loadtxt("results")
        print "The lowest energy in results is: ",
        print data.min()
        print "The index of the lowest energy in results is: ",
        n = data.argmin()
        print data.argmin()
        flag = 1
    if os.path.exists("geo"):
        labels = []
        f = open("geo", "r")
        for i in f:
            if i.startswith("DESCRP"):
                tokens = i.strip().split()
                labels.append(tokens[1])
        if flag == 1:
            print "possibel reference is: ",
            print labels[n]
    print "-----------------------end----------------------"
        



if len(sys.argv) < 3:
    usage()
    preprocess()
else:
    rname = sys.argv[1]
    rener = float(sys.argv[2])
    if len(sys.argv) == 3:
        print "Using default energy unit (kCal/mol)"
        unit = 1
    else:
        if sys.argv[3] == "ev":
            unit = 23
        elif sys.argv[3] == "hartree":
            unit = 623
        elif sys.argv[3] == "kcal":
            unit = 1 
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
        print len(conf), len(ener)
        print "Error!"
    o.close()
    
