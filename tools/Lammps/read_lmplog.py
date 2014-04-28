""" Read lammps log file
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def usage():
    print """python read_lmplog logfile key
    """

def getEnerTerms(logfile):
    keys = []

    f = open(logfile, "r")

    for i in f:
        if "---------------- Step" in i:
            break

    for i in f:
        if "---------------- Step" in i:
            break
        elif "Loop time" in i:
            break
        else:
            tokens = i.strip().split()
            n =  len(tokens)/3
            for i in range(n):
                keys.append(tokens[i*3])
    f.close()

    print "Terms you can select: "
    for i in keys:
        print i

def getData(logfile, key):
    data = []
    step = []
    f = open(logfile, "r")
    for i in f:
        if "---------------- Step" in i:
            break
    for i in f:
        if key in i:
            n = 0
            tokens = i.strip().split()
            for j in tokens:
                if key == j:
                    data.append(float(tokens[n+2]))
                elif "Step" == j:
                    step.append(int(tokens[n+2]))
                n += 1

    assert len(data) > 0, "No terms in log file "
    
    if len(step) == 0:
        x = np.arange(len(data))
    else:
        x = np.array(step)
    ave = np.average(data[len(data)/10:])
    std = np.std(data[len(data)/10:])

    print ave, std

    """
    data = np.array(data)
    textx = (np.average(x) + np.min(x))/2
    texty = (np.average(data) + np.max(data))/2
    plt.text(textx, texty, "%.3f(%.3f)"%(ave, std), size="x-large")
    plt.plot(x, data, "-o")
    plt.plot(x, [ave]*len(x), "--")
    if len(step) > 0:
        plt.xlabel = "Simulation Steps"
    plt.ylabel(key)
    plt.savefig("%s.eps"%key)
    plt.show()
    """

def main():
    if len(sys.argv) < 2:
        usage()
    else:
        logfile = sys.argv[1]
        if len(sys.argv) > 2:
            key = sys.argv[2]
            getData(logfile, key)
        else:
            getEnerTerms(logfile)

main()

