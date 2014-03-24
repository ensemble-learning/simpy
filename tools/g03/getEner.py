#!/usr/bin/env python
"""Get following info: 
   energy, charge, configuration
   from g03 log file. And generate freq calculation input
   automatically.
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

from g03 import G03LogConf, G03tools
from output_conf import toGjf

def main():
    """get qm info
    """
    if len(sys.argv) < 2:
        print "g03.py logfile"
    else:
        for i in sys.argv[1:]:
            geo = i.split(".")[0]
            # get configuration
            a = G03LogConf(i)
            # get charge and energy (HF and Zero)
            b = G03tools(i)
            # parse the configuration to standard format
            c = a.parser()
            hf, zpe = b.getEnergy()
            o = open("qm_ener.dat", "w")
            o.write("%-15s HF %.8f\n"%(geo, hf))
            o.write("%-15s ZPE %.8f\n"%(geo, zpe))
            o.close()
            charges = b.getCharge()
            o = open("qm_charge.dat", "w")
            for j in charges:
                o.write("%-15s%6s%6s%12s  # %s\n"%( geo, "1", j[0], j[2], j[1]))
            o.close()
            # freq calculation if quit normally
            if a.stat == 1:
                method = " ".join([j.strip() for j in c.methods])
                tokens = method.strip().split()
                for n in range(len(tokens)):
                    if "opt" in tokens[n]:
                        tokens[n] = "freq"
                    elif "geo" in tokens[n]:
                        tokens[n] = ''
                c.methods = [" ".join(tokens)]
                c.options = ["%mem=800MB"]
                toGjf(c, "freq.gjf")
            else:
                runfile = os.path.join(os.getcwd(), i)
                print "Warning: unfinished calculation!: %s"%runfile
                method = " ".join([j.strip() for j in c.methods])
                tokens = method.strip().split()
                for n in range(len(tokens)):
                    if "geo" in tokens[n]:
                        tokens[n] = ''
                c.methods = [" ".join(tokens)]
                toGjf(c, "conti.gjf")

if __name__ == "__main__":
    main()
