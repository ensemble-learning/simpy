#!/usr/bin/env python

"""
get energy terms from fort.73
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from plot99 import parse_fort99

class Fort73():
    def __init__(self, ):
        self.e_total = []
        self.e_bond = []
        self.e_atom = []
        self.e_lp = []
        self.e_mol = []
        self.e_val = []
        self.e_coa = []
        self.e_hbo = []
        self.e_tors = []
        self.e_conj = []
        self.e_vdw = []
        self.e_coul = []
        self.e_charge = []
        self.parse_fort73()

    def parse_fort73(self,):
        # read the fort.73 file
        assert os.path.exists("fort.73")

        lines = []
        f = open("fort.73", "r")
        for i in f:
            lines.append(i)
        f.close()

        data = []
        for i in range(len(lines)/3):
            tokens = lines[3*i+2].strip().split()
            data.append([float(j) for j in tokens])

        # parse the energy terms
        ener = np.array(data)

        for i in ener:
            self.e_total.append(np.sum(i[1:]))

        ener_sep = ener.transpose()
        self.e_bond = ener_sep[1]
        self.e_atom = ener_sep[2]
        self.e_lp = ener_sep[3]
        self.e_mol = ener_sep[4]
        self.e_val = ener_sep[5]
        self.e_coa = ener_sep[6]
        self.e_hbo = ener_sep[7]
        self.e_tors = ener_sep[8]
        self.e_conj = ener_sep[9]
        self.e_vdw = ener_sep[10]
        self.e_coul = ener_sep[11]
        self.e_charge = ener_sep[12]

def main():

    # read energy data from fort.73
    ener = Fort73()
    # read the qm data
    assert os.path.exists("fort.99")
    reax, qm = parse_fort99()
    qm = np.array(qm)

    # get bond length if any
    if os.path.exists("bonds"):
        x = np.loadtxt("bonds")
    else:
        print "warning: no bonds were read!!!"
        x = np.arange(len(ener.e_total))

    e_total = ener.e_total - ener.e_total[-1]
    e_bond = ener.e_bond - ener.e_bond[-1]
    e_atom = ener.e_atom - ener.e_atom[-1]
    e_vdw = ener.e_vdw - ener.e_vdw[-1]
    qm = qm - qm[-1]
    plt.plot(x, e_total, '-o', label="Potential energy")
    plt.plot(x, e_bond + e_atom, '-o', label="Bond energy")
    plt.plot(x, e_vdw, '-o', label="vdW energy")
    plt.plot(x, qm, lw=2, color="black", label="QM")

    plt.xlabel(r"Bond length ($\AA$)", size="x-large")
    plt.ylabel("Potential energy (kcal/mol)", size="x-large")
    plt.legend()
    plt.savefig("bondEnergy.eps")

    plt.show()

if __name__ == "__main__":
    main()
