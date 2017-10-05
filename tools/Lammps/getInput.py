""" This code is used to generate reaxFF lammps_input file
according to lammps.data (also generated from simpy
"""
import sys
import os
import socket
import argparse

LIB = ''

print socket.gethostname() 

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "fermion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant3":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant1":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "zwicky":
    LIB = "/home/tcheng/Soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420":
    LIB = "/home/tao/Soft/simpy/lib"
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/simpy/lib"

sys.path.insert(0 , LIB)

from ffield import Ffield

#C H O N
#FF = {"C": 1, "H": 2, "O": 3, "N": 4, "Ca":4, "Al":6}
FF = {}
MASS = {12.011:"C", 14.007: "N", 15.999:"O", 1.0079:"H", 40.078:"Ca",\
        26.982:"Al", 28.086:"Si", 35.453:"Cl", 47.867:"Ti",  6.941:"Li",
        30.974:"P", 32.065:"S", 72.64:"Ge", 10.811:"B", 63.546:"Cu", 
        58.693:"Ni", 195.08:"Pt", 50.942:"V", 95.94:"Mo", 92.906:"Nb",
       127.6:"Te", 22.990:"Na", 69.723:"Ga", 58.933:"Co", 39.948:"Ar", 
        137.327: "Ba", 88.906:"Y", 91.224:"Zr", 126.90:"I",
        107.87: "Ag", 101.07: "Ru", 196.97: "Au", 55.845:"Fe",
        183.84: "W", 65.39: "Zn",
        }

"""
  "Cl":35.453, "P":30.974, "S":32.065, "Ge":72.64, "GE":72.64, "B": 10.811, "ZR":91.224, "V": 50.942, "Mo":95.94,
                "MO": 95.94, "Te": 127.6, "TE": 127.6, "Nb":92.906, "NB":92.906, "Cu":63.546, "Ni": 58.693, "NI": 58.693,
                "Pt": 195.08, "PT":195.08}
"""

def usage():
    print """python genInput.py type
    type: NVT, MIN
    """

def getElements(args):
    assert os.path.exists("ffield")
    if args.lg:
        a = Ffield("ffield", 1 )
    else:
        a = Ffield("ffield", 0 )
    counter = 1
    for i in a.elements:
        FF[i] = counter
        counter += 1

def parseData(fname="lammps.data"):
    """parse the data file to get the mass infor
    """
    m = []
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("Masses"):
            break
    for i in f:
        if len(i.strip()) == 0:
            pass
        elif i.strip().startswith("Atoms"):
            break
        else:
            tokens = i.strip().split()
            amass = float(tokens[1])
            m.append(amass)
    f.close()
    return m

def main(args):
    """ generate the lammps input file
    """
    getElements(args)
    m = parseData()
    ty = []
    elem = []
    for i in m:
        #i = int(i*1000)/1000.0
        at = MASS[i]
        elem.append(at)
        ff = FF[at]
        ty.append("%d"%ff)

    rtype = args.type
    lg = 0
    if args.lg:
        lg = 1
    if rtype == "NVT":
        lines = NVT
    elif rtype == "MIN":
        lines = MIN
    elif rtype == "NPT":
        lines = NPT
    elif rtype == "MIN_CELL":
        lines = MIN_CELL
    elif rtype == "RERUN":
        lines = RERUN 
    elif rtype == "RERUN_REAX":
        lines = RERUN_REAX
    elif rtype == "RERUN_LJ":
        lines = RERUN_LJ

    print "processing %s simulation......"%rtype

    
    if lg:
        lines = lines.replace("%reax_potential%", "reax/c control.reaxc lgvdw yes")
    else:
        lines = lines.replace("%reax_potential%", "reax/c control.reaxc")
    
    if args.lammps2012:
        lines = lines.replace("%ffield_atoms%", " ".join(ty))
    else:
        lines = lines.replace("%ffield_atoms%", " ".join(elem))

    lines = lines.replace("%elements%", " ".join(elem))

    if rtype == "RERUN_LJ":
        lines = RERUN_LJ
        tmp = ""
        for i in range(len(elem)):
            tmp += "pair_coeff     %d     %d      0.0100      1.0000\n"%(i+1, i+1)
        print tmp
        lines = lines.replace("%pair_coeff%", tmp)

    o = open("lammps_input", "w")
    o.write(lines)
    o.close()
    
    o = open("control.reaxc", "w")
    o.write(CONTROL)
    o.close()   
    
    o = open("plumed.dat", "w")
    o.write(PLUMED)
    o.close()   

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("type", default="MIN", nargs="?", help="simulation type")
    parser.add_argument("-lg", action="store_true", help="using lg type ffield")
    parser.add_argument("-lammps2012", action="store_true", help="using lg type ffield")
    parser.add_argument("-template", nargs=1, help="assign lammps input template")
    #parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    #parser.add_argument("-b", nargs=2, type=int, help="get the bond distance between a1, a2, a3")
    #parser.add_argument("-a", nargs=3, type=int,help="get the angle of a1-a2-a3")
    #parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")
    args = parser.parse_args()

    if args.template:
        template_path = args.template[0]
        template_path = os.path.join(template_path, "template")
        sys.path.insert(0 , template_path)
    from template import *

    main(args)


