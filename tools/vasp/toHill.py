#! /usr/bin/env python

"""
Convert the HILLSPOT file (VASP) to HILLS (plumed)

@author: sanli
@author: sanli@llc.com
@version: 1.0
@note: only support the Collective Variables up to 2.
"""

import os

class Hills():
    def __init__(self,):
        """ calss for Hill information
        """
        self.nd = 1
        """
        @ivar:dimensional of collective variables
        """
        self.hills = []
        """
        @ivar:hill information
        """

def read_hills(hills):
    """
    Read Hills from HILLSPOT (VASP)
    """
    f = open("HILLSPOT", "r")
    for i in f:
        tokens = i.strip().split()
        if len(tokens) == 3:
            hills.nd = 1
            hills.hills.append([float(j) for j in tokens])
        elif len(tokens) == 4:
            hills.nd = 2
            hills.hills.append([float(j) for j in tokens])
    f.close()

def output_HILLS_1d(hills):
    """
    Output Hills to plumbed format. 1-D case.
    """
    o = open("HILLS", "w")
    o.write("#! FIELDS time d1 sigma_d1 height biasf\n")
    o.write("#! SET multivariate false\n")
    nstep = 1
    for i in hills.hills:
        rc = i[0]
        sigma = i[2]
        gw = i[1]
        o.write("%8d %12.6f %12.6f %12.6f %4d\n"
                %(nstep, rc, sigma, gw, 1))
        nstep += 1
    o.close()

def output_HILLS_2d(hills):
    """
    Output Hills to plumbed format. 2-D case.
    """
    o = open("HILLS", "w")
    o.write("#! FIELDS time d1 d2 sigma_d1 sigma_d2 height biasf\n")
    o.write("#! SET multivariate false\n")
    nstep = 1
    for i in hills.hills:
        rc1 = i[0]
        rc2 = i[1]
        sigma = i[3]
        gw = i[2]
        o.write("%8d %12.6f %12.6f %12.6f %12.6f %12.6f %4d\n"
                %(nstep, rc1, rc2, sigma, sigma, gw, 1))
        nstep += 1
    o.close()

def run_fe():
    """
    Run free-energy analysis
    """
    os.system("plumed sum_hills --hills HILLS")

def plot_fe():
    import numpy as np
    import matplotlib.pyplot as plt
    
    data = np.loadtxt("fes.dat")
    data = data.transpose()
    plt.plot(data[0], data[1])
    plt.show()

def main():

    hills = Hills()
    log = open("meta.log", "w")

    read_hills(hills)
    if hills.nd == 1:
        output_HILLS_1d(hills)
    elif hills.nd == 2:
        output_HILLS_2d(hills)
    run_fe()
    plot_fe()


if __name__ == "__main__":
    main()
