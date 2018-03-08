"""
@ref: https://plot.ly/python/peak-fitting/
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy
import peakutils

from scipy import signal

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def peak_fit(fname):
    data = np.loadtxt(fname)
    data = data.transpose()

    indices = peakutils.indexes(data[1], thres=0.0, min_dist=0.05)

    peaks = [[],[]]
    for i in range(len(indices)):
        peaks[0].append(data[0][indices[i]])
        peaks[1].append(data[1][indices[i]])

    if 0:
        plt.plot(data[0], data[1])
        plt.scatter(peaks[0], peaks[1]) 
        plt.show()

def main():
    if len(sys.argv) > 1:
        fname = peak_fit(sys.argv[1])
    
if __name__ == "__main__":
    main()
