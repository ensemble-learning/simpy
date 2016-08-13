"""
Get free energy information from VASP slow-growth simulation
Plot the free energy plot derived from MD simulation
"""

import sys
import argparse
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

NSTART = 20

class Params():
    def __init__(self,):
	self.nconst = 1
	self.nxy = 1
	self.nco = 1
	self.nstart_flag = 0
	self.nstart = 0
	self.nend_flag = 0
	self.nend = 0

def read_report(p):
    cvs = []
    dA = []
    fvall = []
    mc = []

    for i in range(p.nconst):
        cvs.append([])
        dA.append([])
        fvall.append([])

    f = open("REPORT", "r")

    nx = 0
    ny = 0
    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("cc> R"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> X"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> Z"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> S"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> C"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> A"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> T"):
            fv1 = float(tokens[2])
            n = nx%p.nconst
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("mc> S"):
            fv1 = float(tokens[2])
            mc.append(fv1)
        elif i.strip().startswith("b_m>"):
            fv1 = float(tokens[1]) # lambda
            fv2 = float(tokens[2]) # |z|^(-1/2)
            fv3 = float(tokens[3]) # GkT
            fv4 = float(tokens[4]) # |z|^(-1/2)*(lambda+GkT)
            n = ny%p.nconst
            fvall[n].append([fv1, fv2, fv3, fv4])
            val = fv4/fv2
            #print val, fv1, fv2, fv3, fv4
            dA[n].append(val)
            ny += 1
    
    for i in range(p.nconst):
        x, p1, p2 = [], [], []
        o = open("fe_data_c%d.dat"%i, "w")
        o.write("#%12s%12s%12s%12s\n"%("colvar", "|z|^(-1/2)", "|z|^(-1/2)*(lambda+GkT)", "dA"))
        for j in range(len(cvs[i])):
            x.append(cvs[i][j])
            p1.append(fvall[i][j][1])
            p2.append(fvall[i][j][3])
            o.write("%12.4f%12.4f%12.4f%12.4f\n"%(cvs[i][j], fvall[i][j][1], fvall[i][j][3], dA[i][j]))
        o.close()

	nstart = NSTART
	nend = len(x)
	if p.nstart_flag == 1:
	    nstart = p.nstart
            if nstart >= nend:
		sys.stderr.write("Start step beyond the simulation step!\n")
		sys.exit()

	if p.nend_flag == 1:
	    if nend < p.nend:
                s.stderr.write("End step beyond the simulation step!\n")
                s.stderr.write("Reset to the simulation end step\n")
	    else:
                nend = p.nend
	
        if nend > nstart:
            p1 = np.array(p1[nstart:nend])
            p2 = np.array(p2[nstart:nend])
            x_ave = np.average(x[nstart:nend])
            y_ave = np.average(p2)/np.average(p1)
            std = np.std(p2)/np.average(p1)
    
        o = open("ave_c%d.dat"%i, "w")
        o.write("%12.4f%12.4f%12.4f\n"%(x_ave, y_ave, std))
        o.close()
    
    return cvs, dA, mc

def plot_data(cvs, dA, mc, p):
    """
    Plot the data
    """
    nxy = p.nxy -1 
    x, y = cvs[nxy], dA[nxy]
    y_int = integrate.cumtrapz(y, x, initial=0)

    o = open("dA.dat", "w")
    for i in range(len(x)):
        o.write("%12.4f%12.4f\n"%(x[i], y[i]))
    o.close()

    o = open("fe.dat", "w")
    for i in range(len(x)):
        o.write("%12.4f%12.4f\n"%(x[i], y_int[i]))
    o.close()
    
    fig = plt.figure(figsize=(8,12))
    ax = fig.add_subplot(2,1,1)
    ax.plot(x, y)
    ax.plot([x[0],x[-1]], [0,0], ls="--")
    ax.set_xlim([x[0], x[-1]])
    ax.set_ylabel(r"$\frac{\partial A}{\partial \xi}$ (ev/A)", size="x-large")
    ax = fig.add_subplot(2,1,2)
    ax.plot(x, y_int, lw=3)
    ax.grid(True)
    ax.set_xlim([x[0], x[-1]])
    ax.set_ylabel("FE (ev)", size="x-large")
    ax.set_xlabel("CV", size="x-large")
    plt.savefig("fe.png")
    plt.show()

def main():
    """
    Extract the information from VASP calculation to free energy
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot", action="store_true", help="plot the data")
    parser.add_argument("-nconst", type=int, nargs=1, help="number of constraints")
    parser.add_argument("-nxy", type=int, nargs=1, help="constraint id to plot")
    parser.add_argument("-nco", type=int, nargs=1, help="variable to observe")
    parser.add_argument("-begin", type=int, nargs=1, help="step to begin with")
    parser.add_argument("-end", type=int, nargs=1, help="step to end up with")
    args = parser.parse_args()

    p = Params()
    if args.nconst:
        p.nconst = args.nconst[0]
    if args.nxy:
        p.nxy = args.nxy[0]
    if args.begin:
        p.nstart_flag = 1
        p.nstart = int(args.begin[0])
    if args.end:
        p.nend_flag = 1
        p.nend = int(args.end[0])
    if args.nco:
        p.nco = int(args.nco[0])

    cvs, dA, mc = read_report(p)

    if args.nco:
        o = open("co.dat", "w")
        for i in range(len(mc)):
            if i%p.nco == 0:
                o.write("\n")
                o.write("%12d"%(int(i/p.nco)))
            o.write("%12.4f"%(mc[i]))
        o.close()

    if args.plot:
        plot_data(cvs, dA, mc, p)

if __name__ == "__main__":
    main()
