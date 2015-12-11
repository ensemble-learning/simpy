"""
Calculate the coordination number according to cutoff functions.
"""

import math
import sys
import ConfigParser

class Params():
    """
    Parameters for coordination calculation.
    """
    def __init__(self, ):
        self.r_ref_0 = 0.0
        self.r_ref_1 = 0.0
        self.r_ref_2 = 0.0
        self.n_atoms = 0
        self.cut_off_func = 1
        self.r_cut_range = 1.0

INP = """[GLOBAL]
r_ref_0 = 2.788
r_ref_1 = 3.944
r_ref_2 = 4.830

n_atoms = 8183
# cut_off_func = 1(cut-off) 2(cn function)
cut_off_func = 1
r_cut_range = 1.5
n_cut_off = 4

"""

def cal_cn(tokens):
    """
    Calculate coordination number using switch function.
    """
    r0 = R0
    d0 = D0
    n = N
    m = M

    a1 = int(tokens[0])
    a2 = int(tokens[1])
    r12 = float(tokens[2])
    p1 = math.pow((r12 - d0)/r0, 6)
    p2 = math.pow((r12 - d0)/r0, 12)
    
    s = (1-p1)/(1-p2)
    return s

def cal_cn_cutoff(tokens, p):
    """
    Calculate coordination number using cutoff.
    """
    rcut = p.r_ref_0 * p.r_cut_range
    #print rcut, p.r_ref_0, p.r_cut_range
    if rcut >= p.r_ref_1:
        rcut = p.r_ref_1
    a1 = int(tokens[0])
    a2 = int(tokens[1])
    r12 = float(tokens[2])
    s = 0
    if r12 < rcut:
        s = 1
    return a1, a2, s

def read_evdw(p):
    """
    Read the vdw file.
    """
    cn = [0] * p.n_atoms
    f = open("lammps.evdw.0", "r")

    counter = 0
    for i in f:
        if counter >= 2:
            tokens = i.strip().split()
            a1, a2, s = cal_cn_cutoff(tokens, p)
            cn[a1-1] += s
            cn[a2-1] += s
        counter +=1
    f.close()
    return cn

def output_ndx(cn, p, log):
    """
    Output the ndx file
    """
    max = 0
    n = 0
    hist = [0]*24
    o = open("low.dat", "w")
    ndx = open("ndx.dat", "w")
    for i in range(len(cn)):
        if cn[i] <= p.n_cut_off:
            n += 1
            o.write("%10d%10.2f\n"%(i+1, cn[i]))
            ndx.write("%8d"%(i+1))
            if n%10 == 0:
                ndx.write(" &\n")
        if cn[i] > max:
            max = cn[i]
        hist[cn[i]-1] += 1
    o.close()
    ndx.close()
    
    log.write("The max coordination number is %d\n"%max)
    log.write("%d atoms with coordination of %d\n"
            %(n, p.n_cut_off))
    log.write("Hist:\n")

    for i in range(max):
        print i
        log.write("%8d%8d\n"%(i+1, hist[i]))

def read_inp(p):
    """Read the input file
    """
    cf = ConfigParser.ConfigParser()
    cf.read("cn.inp")

    s = cf.sections()
    if "GLOBAL" in s:
        o = cf.options("GLOBAL")
        if "r_ref_0" in o:
            tokens = cf.get("GLOBAL", "r_ref_0").strip()
            p.r_ref_0 = float(tokens)
        if "r_ref_1" in o:
            tokens = cf.get("GLOBAL", "r_ref_1").strip()
            p.r_ref_1 = float(tokens)
        if "r_ref_2" in o:
            tokens = cf.get("GLOBAL", "r_ref_2").strip()
            p.r_ref_2 = float(tokens)
        if "n_atoms" in o:
            tokens = cf.get("GLOBAL", "n_atoms").strip()
            p.n_atoms = int(tokens)
        if "cut_off_func" in o:
            tokens = cf.get("GLOBAL", "cut_off_func").strip()
            p.cut_off_func = int(tokens)
        if "r_cut_range" in o:
            tokens = cf.get("GLOBAL", "r_cut_range").strip()
            p.r_cut_range = float(tokens)
        if "n_cut_off" in o:
            tokens = cf.get("GLOBAL", "n_cut_off").strip()
            p.n_cut_off = float(tokens)

def write_inp(p, log):
    log.write("r_ref_0 = %8.4f\n"%params.r_ref_0)
    log.write("r_ref_1 = %8.4f\n"%params.r_ref_1)
    log.write("r_ref_2 = %8.4f\n"%params.r_ref_2)
    log.write("%d of atoms in the system\n"%params.n_atoms)
    log.write("cut_off_func = %d\n"%params.cut_off_func)
    if params.cut_off_func == 1:
        log.write("Using simple cut-off function.\n")
    log.write("cut range = %.4f\n"%(params.r_cut_range))
    log.write("    cut-off = %.4f\n"%(params.r_ref_0 * params.r_cut_range))

def main(p, log):
    # read vdw input
    cn = read_evdw(p)
    # output ndx 
    output_ndx(cn, p, log)

if __name__ == "__main__":
    log = open("cn.log", "w")

    params = Params()
    read_inp(params)
    write_inp(params, log)
    main(params, log)

    log.close()
