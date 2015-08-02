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
        elif i.strip().startswith("b_m>"):
            fv1 = float(tokens[1]) # lambda
            fv2 = float(tokens[2]) # |z|^(-1/2)
            fv3 = float(tokens[3]) # GkT
            fv4 = float(tokens[4]) # |z|^(-1/2)*(lambda+GkT)
            fvall.append([fv1, fv2, fv3, fv4])
            val = fv4/fv2
            dA.append(val)
    
    p1 = 0.0
    p2 = 0.0
    x = 0.0
    n = 0
    for i in range(len(cvs)):
        x += cvs[i]
        p1 += fvall[i][1]
        p2 += fvall[i][3]
        n += 1
    x = x/n
    p1 = p1/n
    p2 = p2/n
    y = p2/p1
    o = open("ave.dat", "w")
    o.write("%12.4f%12.4f\n"%(x, y))
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
