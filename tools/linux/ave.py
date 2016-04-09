#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) > 2:
        fname = sys.argv[1]
        nr = sys.argv[2]

        data = np.loadtxt(fname)
        data = data.transpose()
 
    	f = 0.1
    	nstart = int(f*len(data))

    	plt.plot(data[0], data[1])
    	ave = np.average(data[1][nstart:])
    	print ave
    	plt.show()
    

if __name__ == "__main__":
    main()
