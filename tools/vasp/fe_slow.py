"""
Get free energy information from VASP slow-growth simulation
Plot the free energy plot derived from MD simulation
"""

import argparse
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def read_report():
    cvs = []
    dA = []
    fvall = []
    f = open("REPORT", "r")

    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("cc> R"):
            fv1 = float(tokens[2])
            cvs.append(fv1)
        elif i.strip().startswith("cc> C"):
            fv1 = float(tokens[2])
            cvs.append(fv1)
        elif i.strip().startswith("cc> X"):
            fv1 = float(tokens[2])
            cvs.append(fv1)
        elif i.strip().startswith("cc> S"):
            fv1 = float(tokens[2])
            cvs.append(fv1)
        elif i.strip().startswith("cc> Z"):
            fv1 = float(tokens[2])
            cvs.append(fv1)
        elif i.strip().startswith("b_m>"):
            fv1 = float(tokens[1]) # lambda
            fv2 = float(tokens[2]) # |z|^(-1/2)
            fv3 = float(tokens[3]) # GkT
            fv4 = float(tokens[4]) # |z|^(-1/2)*(lambda+GkT)
            fvall.append([fv1, fv2, fv3, fv4])
            val = fv4/fv2
            dA.append(val)
    
    x, p1, p2 = [], [], []
    o = open("fe_data.dat", "w")
    for i in range(len(cvs)):
        x.append(cvs[i])
        p1.append(fvall[i][1])
        p2.append(fvall[i][3])
        o.write("%12.4f%12.4f%12.4f\n"%(cvs[i], fvall[i][1], fvall[i][3]))
    o.close()

    x = np.array(x)
    if len(x) >= 10:
        nstart = 10
        p1 = np.array(p1[nstart:])
        p2 = np.array(p2[nstart:])
        x = np.average(x[nstart:])
        y = np.average(p2)/np.average(p1)
        std = np.std(p2)/np.average(p1)
    
    o = open("ave.dat", "w")
    o.write("%12.4f%12.4f%12.4f\n"%(x, y, std))
    o.close()
    
    return cvs, dA

def plot_data(cvs, dA):
    x, y = cvs, dA
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
    ax.plot(x, y_int)
    ax.grid(True)
    ax.set_xlim([x[0], x[-1]])
    ax.set_ylabel("FE (ev)", size="x-large")
    ax.set_xlabel("CV", size="x-large")
    plt.savefig("fe.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot", action="store_true", help="plot the data")
    args = parser.parse_args()

    cvs, dA = read_report()
    if args.plot:
        plot_data(cvs, dA)

if __name__ == "__main__":
    main()
