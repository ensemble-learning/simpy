""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import sys, socket
import copy
import os
import math

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
from output_conf import toPdb, toReaxLammps

DELTA_Z1 = 15.0
DELTA_Z2 = 5.0

def write_lammps(ca,c1,c2,a_ref):
    """
    Write lammps constraint file (constaint.inc).
    @param ca: all atoms
    @param c1: fixed atoms
    @param c2: relaxed atoms
    @param a_ref: reference atom
    """
    n1 = len(c1)
    n2 = len(ca)
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

def write_cn_inp(ca, c1, c2, a_ref):
    """
    Write coordination input file (cn.inp).
    @param ca: all atoms
    @param c1: fixed atoms
    @param c2: relaxed atoms
    @param a_ref: reference atom
    """
    o = open("cn.inp", "w")
    o.write("""[GLOBAL]
r_ref_0 = 2.788
r_ref_1 = 3.944
r_ref_2 = 4.830

n_atoms = %d
# cut_off_func = 1(cut-off) 2(cn function)
cut_off_func = 1
r_cut_range = 1.2
n_cut_off = 4
ignore_atoms = %d
"""%(len(ca), len(c1)))
    o.close()

def get_slice(b, id0):
    logfile = "build.log"
    o = open(logfile, "w")
    
    a0 = b.atoms[id0-1]
    id0 = a0.an
    x0, y0, z0 = a0.x

    o.write("Reference atoms is %d.\n"%id0)
    o.write("The coordination of reference atom is %.4f %.4f %.4f.\n"%
            (x0, y0, z0))
    o.write("Cut from z direction.\n")

    zl = b.pbc[2]

    z0 = z0 - int(z0/zl)*zl
    dz1 = DELTA_Z1
    dz2 = DELTA_Z2
    z1 = z0 - dz1
    z2 = z0 + dz1
    z1 = z1 - math.floor(z1/zl)*zl
    z2 = z2 - math.floor(z2/zl)*zl

    b1 = copy.deepcopy(b)

    o.write("Starting from %.4f to %.4f.\n"%(z1, z2))
    o.write("Cut radius is %.4f.\n"%dz1)
    #o.write("Buffer distance is %.4f.\n"%dz2)
    #o.write("Fixed L1 from %.4f to %.4f\n"%(z1, z2))
    #o.write("Fixed L2 from %.4f to %.4f\n"%(z3, z4))
    
    ca = []
    a_ref = []
    
    for i in range(len(b.atoms)):
        id = b.atoms[i].an
        z = b.atoms[i].x[2]
        z = z - math.floor(z/zl)*zl #pbc
        b1.atoms[i].x[2] = z
        if id == id0:
            a_ref.append(b1.atoms[i])
        else:
            if z1 < z2:
                if z >= z1 and z < z2:
                    ca.append(b1.atoms[i])
            else:
                if z >= 0 and z < z2:
                    ca.append(b1.atoms[i])
                elif z >= z1 and z < zl:
                    ca.append(b1.atoms[i])

    # Shift the atoms according to reference atom
    for i in range(len(ca)):
        z = ca[i].x[2] 
        z = z - z1
        z = z - math.floor(z/zl)*zl
        ca[i].x[2] = z

    z = a_ref[0].x[2] 
    z = z - z1
    z = z - math.floor(z/zl)*zl
    a_ref[0].x[2] = z
    o.write("Sort the atom index according to z.\n")
    ca.sort(key=lambda x: x.x[2])

    c1 = []
    c2 = []
    for i in ca:
        z = i.x[2]
        if z <= dz2:
            c1.append(i)
        elif z > 2*dz1 - dz2:
            c1.append(i)
        else:
            c2.append(i)
            
    o.write("Fixed atoms is %d.\n"%(len(c1)))

    ca = ca + a_ref
    o.write("Total number of atoms is %d.\n"%len(ca))

    o.write("Index is as following:\n")
    counter = 0
    for i in ca:
        if counter % 10 == 0:
            o.write("\n")
        o.write("%8d"%i.an)
        counter += 1
    o.write("\n")
    
    b1.atoms = c1 + c2 + a_ref
    b1.pbc[2] = 2*dz1
    
    o.write("Output to pdb.\n")
    toPdb(b1, "input.pdb")
    toReaxLammps(b1)
    o.close()

    return ca, c1, c2, a_ref

def read_data():
    datafile = "lammps.data"
    a = ReaxData(datafile)
    b = a.parser()
    b.assignAtomTypes()
    return b
    
def read_ndx_coords(log):
    """
    Read index and coordiantions from ndx coords file (coords.dat),
    and check overlap.
    @param log: log file
    @return: ndx list
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
                if z_now - z_prev > 2*DELTA_Z1:
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
    log.write("Leaving the overlapped atoms to nex round:\n")
    for i in range(len(res)):
        if i%10 == 0:
            log.write("\n")
        log.write("%8d"%res[i])
    
    return ndx

def main(log):
    """
    Slice driver
    """
    ndx = read_ndx_coords(log)
    data = read_data()

    for i in range(len(ndx)):
        folder = "slice_%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        id0 = ndx[i]
        ca, c1, c2, a_ref = get_slice(data, id0)
        write_lammps(ca,c1,c2,a_ref)
        write_cn_inp(ca,c1,c2,a_ref)
        os.chdir("..")

if __name__ == "__main__":
    log = open("slice.log", "w")
    main(log)
    log.close()
