""" Give a set of esimated bond order parameters
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import poly1d, polyfit

def estimate_bo_params(p1=-0.05):
    data = np.loadtxt("bo.dat")
    data = data.transpose()
    r = data[0]
    bo = data[1]

    x = np.log(r)
    y = np.log(np.log(bo)/p1)

    fit = polyfit(x, y, 1)
    fit_fn = poly1d(fit)

    p2 = fit[0]
    p1 = p1
    r0 = math.exp(-fit[1]/p2)

    print("r0 = ", r0)
    print("pbo1 = ", p1)
    print("pbo2 = ", p2)

    plt.plot(x, y, x,fit_fn(x), '--k')
    plt.show()

def main():
    pass

if __name__ == "__main__":
    if len(sys.argv) > 1:
        p1 = float(sys.argv[1])
        estimate_bo_params(p1)
    else:
        estimate_bo_params()
