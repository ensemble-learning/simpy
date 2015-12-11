""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import sys, socket
import copy
import os

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
elif socket.gethostname() == "giant3":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant1":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif "node" in socket.gethostname():
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from data import ReaxData
from output_conf import toPdb

DELTA_Z1 = 15.0
DELTA_Z2 = 5.0

def get_slice(b, id0):
    logfile = "build.log"
    o = open(logfile, "w")
    
    a0 = b.atoms[id0-1]
    id0 = a0.an
    x0, y0, z0 = a0.x

    dz1 = DELTA_Z1
    dz2 = DELTA_Z2
    z1 = z0 - dz1
    z2 = z1 + dz2
    z3 = z0 + (dz1 -dz2)
    z4 = z0 + dz1
    
    o.write("Reference atoms is %d.\n"%id0)
    o.write("Th coordination of reference atom is %.4f %.4f %.4f.\n"%
            (x0, y0, z0))
    o.write("Cut from z direction.\n")
    o.write("Starting from %.4f.\n"%z0)
    o.write("Cut radius is %.4f.\n"%dz1)
    o.write("Buffer distance is %.4f.\n"%dz2)
    o.write("Fixed L1 from %.4f to %.4f\n"%(z1, z2))
    o.write("Fixed L2 from %.4f to %.4f\n"%(z3, z4))
    
    c1 = []
    c2 = []
    c3 = []
    c4 = []
    
    for i in b.atoms:
        id = i.an
        if id == id0:
            i.x[2] = i.x[2] - z1
            c4.append(i)
        else:
            z = i.x[2]
            if z >= z1 and z < z2:
                i.x[2] = i.x[2] - z1
                c1.append(i)
            elif z >= z2 and z < z3:
                i.x[2] = i.x[2] - z1
                c3.append(i)
            elif z >= z3 and z < z4:
                i.x[2] = i.x[2] - z1
                c2.append(i)
    
    # sort 
    o.write("Sort the atoms according to z.\n")
    c1.sort(key=lambda x: x.x[2])
    c2.sort(key=lambda x: x.x[2])
    c3.sort(key=lambda x: x.x[2])
    
    c = c1 + c2 + c3 + c4
    o.write("Total number of atoms is %d.\n"%(len(c)))
    o.write("Fixed atoms is %d.\n"%(len(c1+c2)))
    o.write("Index is as following:\n")
    
    counter = 0
    for i in c:
        if counter % 10 == 0:
            o.write("\n")
        o.write("%8d"%i.an)
        counter += 1
    o.write("\n")
    
    b1 = copy.deepcopy(b)
    b1.atoms = c
    b1.pbc[2] = 2*dz1
    
    o.write("Output to pdb.\n")
    toPdb(b1, "input.pdb")
    o.close()

    return c1, c2, c3, c4

def read_data():
    datafile = "lammps.data"
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b
    
def read_ndx_coords(log):
    """
    """

    ndx = []
    res = []
    ndxcoordsfile = "coords.dat"
    log.write("Reading %s.\n"%ndxcoordsfile)
    f = open(ndxcoordsfile, "r")
    n = 0
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split()
            id = int(tokens[0])
            z = float(tokens[-1])
            z_now = z
            if n == 0:
                z_prev = z_now
                ndx.append(id)
                n += 1
            else:
                if z_now - z_prev > DELTA_Z1:
                    ndx.append(id)
                    z_prev = z
                else:
                    res.append(id)
    log.write("Dealing with the following atoms:\n")
    for i in range(len(ndx)):
        if i%10 == 0:
            log.write("\n")
        log.write("%8d"%ndx[i])
    log.write("\n")
    log.write("Leaving the overlapped atoms:\n")
    for i in range(len(res)):
        if i%10 == 0:
            log.write("\n")
        log.write("%8d"%res[i])
    
    return ndx

def write_lammps(c1,c2,c3,c4):
    n1 = len(c1+c2)
    n2 = len(c1+c2+c3+c4)
    o = open("constaint.inc", "w") 
    o.write("""group          freeze id <= %d
group          sim1 id %d
group          sim subtract all freeze
fix            901 freeze setforce 0.0 0.0 0.0
"""%(n1, n2))
    o.close()
    o = open("colvar.inp", "w")
    o.write("""colvar {
  name cn
  lowerWallConstant 1000
  lowerWall 2.7
  coordNum {
    group1 { atomNumbers %d}
    group2 { atomNumbersRange 1-%d}
    cutoff 3.288
    expNumer 16
    expDenom 32
  }
}

metadynamics {
  name meta
  colvars cn
  useGrids off
  hillWeight 2.00
  hillWidth 0.2
  newHillFrequency 80
}
"""%(n2, n2-1))
    o.close()

def main(log):
    ndx = read_ndx_coords(log)
    data = read_data()

    for i in range(len(ndx)):
        folder = "slice_%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        id0 = ndx[i]
        c1,c2, c3, c4 = get_slice(data, id0)
        write_lammps(c1,c2,c3,c4)
        os.chdir("..")

if __name__ == "__main__":
    log = open("slice.log", "w")
    main(log)
    log.close()
