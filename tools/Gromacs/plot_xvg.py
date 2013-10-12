import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print "plot_xvg.py [xvgfile]"
    else:
        msdfile = sys.argv[1]
        plot(msdfile)

def plot(msdfile):
    data = []
    f = open(msdfile, "r")
    for i in f:
        if i.strip().startswith("#") or i.strip().startswith("@"):
            if "title" in i:
                title = i[10:].strip()
            pass
        else:
            tokens = i.strip().split()
            if len(tokens) == 0:
                pass
            else:
                data.append([float(j) for j in tokens])
    f.close()
    data = np.array(data)
    data = data.transpose()
    for i in range(1, len(data)):
        plt.plot(data[0], data[i], "-o")
    plt.show()

if __name__ == "__main__":
    main()
                
