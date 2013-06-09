"""
calculate the heat of formation according to the geo file
@note: the data base for standard heat of formation should be updated.
"""

import sys
import socket
LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simupy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib" 

sys.path.insert(0 , LIB)
import os
from geo import Geo

# ref data from ReaxFF 2001
#STD = {"H": 54.3, "C": 218.6 }
STD = {"H": 54.3, "C": 199.6 }


def cal_heatf():
    """calculate the heat of formation as following:
    \delta H_f = E_system + 4RT + POP + n_CI_C + n_HI_H
    POP: contribution of high energy conformations. (default -0.2 kcal/mol)
    I_C: heat increament for C (see STD)
    I_H: heat increament for H (see STD)
    @ref: reaxFF 2001
    """
    # get E_system
    assert os.path.exists("fort.74")
    f =  open("fort.74", "r")
    ener = float(f.readline()[27:37])
    f.close()

    assert os.path.exists("geo")
    fname = "geo"
    a = Geo(fname)
    b = a.parser()
    # get atom map
    b.assignAtomTypes()
    # get element type
    b.assignEleTypes()
    ht = {}

    for i in b.map:
        ht[i[1]] = 0

    for i in b.atoms:
        ener = ener + STD[i.element]

    ener = ener + 4*300*8.314/4184 + (-0.21)

    print b.name, ener

if __name__ == "__main__":
    cal_heatf()
