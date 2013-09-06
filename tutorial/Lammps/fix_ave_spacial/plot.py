import numpy as np
import matplotlib.pyplot as plt

def parse_fix(fname):
    """ parse lammps output file (generated from fix)
    """
    data = []
    n = 4
    for i in range(n):
        data.append([])
    f = open(fname, "r")
    counter = 0
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            if counter == 0:
                pass
            else:
                tokens = i.strip().split()
                if len(tokens) == n:
                    for j in range(n):
                        data[j].append(float(tokens[j]))
            counter += 1
    data = np.array(data)
    return data

def plot_charge():
    data = parse_fix("charge_electrode.dist")
    plt.plot(data[1], data[3]*data[2], lw=3, ls="--", label="Electrode")

    data = parse_fix("charge_electrolyte.dist")
    plt.plot(data[1], data[3]*data[2], lw=3, ls="--", label="Electrolyte")

    data = parse_fix("charge_all.dist")
    plt.plot(data[1], data[3]*data[2], lw=3, label="Total")
    plt.xlabel("Z ($\AA$)", size="x-large")
    plt.ylabel("Charge (e)", size="x-large")
    plt.xlim([0,50])
    plt.legend()
    plt.show()

def plot_number():
    data = parse_fix("number_electrode.dist")
    plt.plot(data[1], data[2], lw=3, ls="--", label="Electrode")

    data = parse_fix("number_li.dist")
    plt.plot(data[1], data[2], lw=3, ls="--", label="Li+")

    data = parse_fix("number_cl.dist")
    plt.plot(data[1], data[2], lw=3, label="Cl-")

    data = parse_fix("number_dmc.dist")
    plt.plot(data[1], data[2], lw=3, label="DMC")
    plt.xlabel("Z ($\AA$)", size="x-large")
    plt.ylabel("Number", size="x-large")
    plt.xlim([0,50])
    plt.legend(loc=2)
    plt.show()

plot_charge()
