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

    sys.stdout.write("Terms you can select: ")
    n = 0
    for i in keys:
        if n%8 == 0:
            sys.stdout.write("\n")
        sys.stdout.write("%-8s\t"%i)
        n += 1
    sys.stdout.write("\n")
            
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

    sys.stdout.write("%8s = %15.2f (%.2f)\n"%(key, ave, std))

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
    plt.savefig("%s.png"%key)
    #plt.show()
    o = open("%s_lammps.csv"%key, "w")
    for i in range(len(x)):
        o.write("%.4f,%.4f\n"%(x[i], data[i]))
    o.close()
"""

def main():
    if len(sys.argv) < 2:
        usage()
    else:
        logfile = sys.argv[1]
        if len(sys.argv) > 2:
            keys = sys.argv[2:]
            for key in keys:
                getData(logfile, key)
        else:
            getEnerTerms(logfile)

main()

