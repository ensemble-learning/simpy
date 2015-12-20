#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    fname = sys.argv[1]

    data = np.loadtxt(fname)
    f = 0.1
    nstart = int(f*len(data))

    plt.plot(data)
    ave = np.average(data[nstart:])
    print ave
    plt.show()
    

if __name__ == "__main__":
    main()
