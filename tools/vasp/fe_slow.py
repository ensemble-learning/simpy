"""
Get free energy information from VASP slow-growth simulation
Plot the free energy plot derived from MD simulation
"""

import argparse
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def read_report(n_const):
    cvs = []
    dA = []
    fvall = []

    print n_const
    for i in range(n_const):
        print i
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
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> X"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> Z"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> S"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> C"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> A"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("cc> T"):
            fv1 = float(tokens[2])
            n = nx%n_const
            cvs[n].append(fv1)
            nx += 1
        elif i.strip().startswith("b_m>"):
            fv1 = float(tokens[1]) # lambda
            fv2 = float(tokens[2]) # |z|^(-1/2)
            fv3 = float(tokens[3]) # GkT
            fv4 = float(tokens[4]) # |z|^(-1/2)*(lambda+GkT)
            n = ny%n_const
            fvall[n].append([fv1, fv2, fv3, fv4])
            val = fv4/fv2
            dA[n].append(val)
            ny += 1
    
    for i in range(n_const):
        x, p1, p2 = [], [], []
        o = open("fe_data_c%d.dat"%i, "w")
        for j in range(len(cvs[i])):
            x.append(cvs[i][j])
            p1.append(fvall[i][j][1])
            p2.append(fvall[i][j][3])
            o.write("%12.4f%12.4f%12.4f%12.4f\n"%(cvs[i][j], fvall[i][j][1], fvall[i][j][3], dA[i][j]))
        o.close()

        if len(x) >= 10:
            nstart = 10
            p1 = np.array(p1[nstart:])
            p2 = np.array(p2[nstart:])
            x_ave = np.average(x[nstart:])
            y_ave = np.average(p2)/np.average(p1)
            std = np.std(p2)/np.average(p1)
    
        o = open("ave_c%d.dat"%i, "w")
        o.write("%12.4f%12.4f%12.4f\n"%(x_ave, y_ave, std))
        o.close()
    
    return cvs, dA

def plot_data(cvs, dA, nxy):
    """
    Plot the data
    """
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
    main code
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot", action="store_true", help="plot the data")
    parser.add_argument("-nconst", type=int, nargs=1, help="number of constraints")
    parser.add_argument("-nxy", type=int, nargs=1, help="constraint id to plot")
    args = parser.parse_args()

    nconst = 1
    if args.nconst:
        nconst = args.nconst[0]
    nxy = nconst - 1
    if args.nxy:
        nxy = args.nxy[0]

    cvs, dA = read_report(nconst)

    if args.plot:
        plot_data(cvs, dA, nxy)

if __name__ == "__main__":
    main()
