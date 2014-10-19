"""
generate Atom.dic and Cutoff.dic
@note: need the ffield.
@todo: need argument parser
"""

import sys
import os
import socket

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from ffield import Ffield
from cons import ELEMENT2MASS

ff = Ffield()
elements = ff.elements

if 1:
    for i in range(len(elements)):
        if len(elements[i]) > 1:
            elements[i] = elements[i][0] + elements[i][1].upper()

# write Atom.dic
o = open("Atom.dic", "w")
for i in elements:
    o.write("%s    %.4f\n"%(i, ELEMENT2MASS[i]))
o.close()

# write Cutoff.dic
o = open("Cutoff.dic", "w")
for i in range(len(elements)):
    for j in range(i, len(elements)):
        o.write("%s    %s    0.3\n"%(elements[i], elements[j]))
o.close()

