import sys
import numpy as np
from numpy import exp, loadtxt, pi, sqrt
import matplotlib.pyplot as plt
from lmfit import Model

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (sqrt(2*pi) * wid)) * exp(-(x-cen)**2 / (2*wid**2))


def gaussian_fit(fname):
    fname = sys.argv[1]

    data = np.loadtxt(fname)
    data = data.transpose()
    x = data[0]
    y = data[1]

    gmodel = Model(gaussian)
    result = gmodel.fit(y, x=x, amp=120, cen=9.1, wid=0.2)
    print(result.fit_report())
    #plt.plot(data[0], data[1])
    #plt.plot(x, result.init_fit, 'k--')
    #plt.plot(x, result.best_fit, 'r-')
    #plt.show()

def main():
    if len(sys.argv) > 1:
        fname = sys.argv[1]
        gaussian_fit(fname)
    
if __name__ == "__main__":
    main()

