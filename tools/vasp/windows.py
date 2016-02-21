"""
Get windows from meta dynamics simulations
for TI calculations
"""
    
import os
import sys
import numpy as np
import ConfigParser

INP = """[PARAMS]
TOSTEP = 1000
x0 = 1.4
x1 = 3.5
n_bins = 15
"""

class Params():
    def __init__(self,):
        self.TOSTEP = 1
        self.x0 = 1.0
        self.x1 = 5.0
        self.n_bins = 15
        self.dx = (self.x1 - self.x0)/self.n_bins

def write_inp():
    o = open("test.inp.0", "w")
    o.write(INP)
    o.close()
        
def read_inp(p):
    """Read the input file
    """
    cf = ConfigParser.ConfigParser()
    cf.read("windows.inp")

    s = cf.sections()
    if "PARAMS" in s:
        o = cf.options("PARAMS")
        if "TOSTEP" in o:
            tokens = cf.get("PARAMS", "TOSTEP").strip()
            p.tostep = int(tokens)
        if "x0" in o:
            tokens = cf.get("PARAMS", "x0").strip()
            p.x0 = float(tokens)
        if "x1" in o:
            tokens = cf.get("PARAMS", "x1").strip()
            p.x1 = float(tokens)
        if "n_bins" in o:
            tokens = cf.get("PARAMS", "n_bins").strip()
            p.n_bins = int(tokens)
    
    p.dx = (p.x1 - p.x0)/p.n_bins
    
    sys.stdout.write("x0 = %.4f\n"%(p.x0))
    sys.stdout.write("x1 = %.4f\n"%(p.x1))
    sys.stdout.write("n_bins = %d\n"%(p.n_bins))
    sys.stdout.write("dx = %.4f\n"%(p.dx))

def read_colvars():
    f = open("COLVARS", "r")
    data = [[],[]]
    for i in f:
        if i.strip().startswith("PLUM"):
            pass
        elif i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split()
            data[0].append(float(tokens[0]))
            data[1].append(float(tokens[1]))
    
    return data

def write_bins(p, data):
    colvar_values = []
    colvar_n = []
    others_values = []
    others_n = []
    for i in range(p.n_bins):
        colvar_values.append([])
        colvar_n.append([])
    
    for i in range(len(data[0])):
        colvar = data[1][i]
        n = int((colvar - p.x0)/p.dx)
        an = int(data[0][i] * 1000)
        print n
        if n >= 0 and n < p.n_bins:
            colvar_values[n].append(colvar)
            colvar_n[n].append(an)
        else:
            others_values.append(colvar)
            others_n.append(an)
        
    o = open("windows.dat", "w")
    for i in range(p.n_bins):
        o.write("%d %.4f\n"%(colvar_n[i][0], colvar_values[i][0]))
    o.close()

def main():
    write_inp()
    p = Params()
    if os.path.exists("windows.inp"):
        read_inp(p)
        data = read_colvars()
        write_bins(p, data)

if __name__ == "__main__":
    main()
