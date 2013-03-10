import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def parseLog(fname):
    """parse the Lammps Log file
    """
    data = []
    keys = []
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("Step"):
            keys = i.strip().split()
            print keys
            break
    for i in f:
        if i.strip().startswith("Loop time"):
            break
        else:
            tokens = i.strip().split()
            data.append([float(j) for j in tokens])
    return keys, data

def main(keys, data):
    data = np.array(data)
    data = data.transpose()
    print "Please select the term:"
    for i in range(1, len(keys)):
        print "%-6d"%i, keys[i]
    print "your choice is:",
    t = int(raw_input())
    print "Average %s = %.4f +/- %.4f"%(keys[t], np.average(data[t]),\
                                        np.std(data[t]))
    # output the data
    o = open("%s.xvg"%keys[t], "w")
    for i in range(len(data[0])):
        o.write("%-12d%-16.4f\n"%(data[0][i], data[t][i]))
    o.close()

    # plot the data
    print "Plot the data?(y):",
    if (raw_input() == "y"):
        plt.plot(data[0], data[t])
        plt.show()

def analysis(fname, keyword):
    data = []
    f = open(fname, "r")
    o = open("%s.xvg"%keyword, "w")
    for i in f:
        if keyword in i:
            tokens = i.strip().split("=")
            for j in range(len(tokens)):
                if keyword in tokens[j]:
                    a = tokens[j+1].strip().split()[0]
                    o.write("%s\n"%a)
                    data.append(float(a))
    o.close()
    f.close()
    data = np.array(data)
    start = int(len(data) * 0.1)
    print np.average(data[start:])

if __name__ == "__main__":
    
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
        keys, data = parseLog(fname)
        main(keys, data)
    else:
        if os.path.exists("log.lammps"):
            fname = "log.lammps"
            print "Using default log file log.lammps"
            keys, data = parseLog(fname)
            main(keys, data)
        else:
            print "readLog.py logfile"
    
    
    
