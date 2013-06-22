#!/usr/bin/env python

"""
@ref: J. Phys. Chem. A 2001, 105, 9396-9409
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt

class FfieldEq():
    def __init__(self, ffield="output.ff"):
        self.eq2 = []
        self.eq4 = []
        self.eq6 = []
        self.eq23 = []
        self.eq26 = []
        self.read(ffield)
    def read(self,ffield):
        assert os.path.exists(ffield)
        f = open(ffield, "r")
        for i in f:
            if i.strip().startswith(r"#"):
                break
        for i in f:
            if i.strip().startswith(r"#"):
                break
            self.eq2.append(i.strip().split())
        for i in f:
            if i.strip().startswith(r"#"):
                break
            self.eq4.append(i.strip().split())
        for i in f:
            if i.strip().startswith(r"#"):
                break
            self.eq6.append(i.strip().split())
        for i in f:
            if i.strip().startswith(r"#"):
                break
            self.eq23.append(i.strip().split())
        for i in f:
            if i.strip().startswith(r"#"):
                break
            self.eq26.append(i.strip().split())

class Bond():
    def __init__(self,):
        self.a1 = '' #atom1
        self.a2 = '' #atom2
        self.p_bo1 = 0.0
        self.p_bo2 = 0.0
        self.p_bo3 = 0.0
        self.p_bo4 = 0.0
        self.p_bo5 = 0.0
        self.p_bo6 = 0.0
        self.r_o = 0.0
        self.r_pi = 0.0
        self.r_pipi = 0.0

        self.lamda_1 = 0.0
        self.lamda_2 = 0.0
        self.lamda_3 = 0.0
        self.lamda_4 = 0.0
        self.lamda_5 = 0.0

        self.vali = 0.0
        self.valj = 0.0

        self.D_e = 0.0
        self.p_be1 = 0.0
        self.p_be2 = 0.0
        self.D_ep = 0.0
        self.D_epp = 0.0
        # vdW parameters
        self.rcut = 1
        self.D_ij = 0.0
        self.alpha_ij = 0.0
        self.r_vdw = 1.0
        self.p_vdw1 = 1.0
        self.gamma_w = 1.0
        # charge parameters
        self.rcut = 4.5
        self.q1 = 1
        self.q2 = -1
        self.gamma_rij = 1.0
        # inner wall parameters
        self.D_inner = 0.0
        self.alpha_inner = 0.0
        self.r_inner = 0.01

        self.r = 0.0
        self.bo = 0.0
        self.bo_s = 0.0
        self.bo_p = 0.0
        self.bo_pp = 0.0
        self.be = 0.0
        self.be_s = 0.0
        self.be_p = 0.0
        self.be_pp = 0.0

    def bond_order(self, r):
        """
        Calculate the uncorrected bond order.
        """
        #print self.r_o, self.p_bo2, self.p_bo1
        self.r = r

        if self.r_o < 0.01:
            print "Warning: too small r_o! Set r = 0.01"
            self.r_o = 0.01

        bo = 0
        bpow = math.pow((r/self.r_o), self.p_bo2)
        bexp = math.exp(self.p_bo1*bpow)
        bo += bexp
        self.bo_s = bexp

        # maybe not safe enouth
        if self.p_bo3 != 0.0 and self.r_pi > 0:
            bpow2 = math.pow((r/self.r_pi), self.p_bo4)
            bexp2 = math.exp(self.p_bo3*bpow2)
            bo += bexp2
            self.bo_p = bexp2
        if self.p_bo5 != 0.0 and self.r_pipi > 0:
            bpow3 = math.pow((r/self.r_pipi), self.p_bo6)
            bexp3 = math.exp(self.p_bo5*bpow3)
            bo += bexp3
            self.bo_pp = bexp3
        self.bo = bo

    def bo_corr(self, bop, deltai, deltaj):
        sum = -self.lamda_3*(self.lamda_4*bop*bop - deltai)
        bexp = math.exp(sum + self.lamda_5) 
        f4 = 1/(1 + bexp)

        sum2 = -self.lamda_3*(self.lamda_4*bop*bop - deltaj)
        bexp2 = math.exp(sum2 + self.lamda_5) 
        f5 = 1/(1 + bexp2)
        
        expi = math.exp(-self.lamda_1*deltai)
        expj = math.exp(-self.lamda_1*deltaj)
        f2 = expi + expj

        expi2 = math.exp(-self.lamda_2*deltai)
        expj2 = math.exp(-self.lamda_2*deltaj)
        f3 = -1/self.lamda_2*math.log(0.5*(expi2 + expj2))
        
        f1 = (self.vali+f2)/(self.vali+f2+f3)
        f1 += (self.valj+f2)/(self.valj+f2+f3)
        f1 = f1/2.0
        return f1, f2, f3, f4, f5

    def bond_energy(self):
        """
        Calculate Bond Order
        """
        print self.D_e
        bo_s = self.bo_s
        bo_p = self.bo_p
        bo_pp = self.bo_pp
        ebond = 0
        bpow = 1 - math.pow(bo_s, self.p_be2)
        bexp = math.exp(self.p_be1*bpow)
        ebond = -self.D_e * bo_s * bexp
        self.be_s = ebond
        if self.D_ep > 0:
            ebond += -self.D_ep*bo_p
            self.be_p = ebond
        if self.D_epp > 0:
            ebond += -self.D_epp*bo_pp
            self.be_pp = ebond
        self.be = ebond

    def tap(self, r):
        t = 0
        tap7 = 20.0/math.pow(self.rcut, 7)
        tap6 = -70.0/math.pow(self.rcut, 6)
        tap5 = 84.0/math.pow(self.rcut, 5)
        tap4 = -35.0/math.pow(self.rcut, 4)
        tap3 = 0.0
        tap2 = 0.0
        tap1 = 0.0
        tap0 = 1.0
        t += tap7*math.pow(r, 7)
        t += tap6*math.pow(r, 6)
        t += tap5*math.pow(r, 5)
        t += tap4*math.pow(r, 4)
        t += tap3*math.pow(r, 3)
        t += tap2*math.pow(r, 2)
        t += tap1*math.pow(r, 1)
        t += tap0
        return t

    def vdw(self, r):
        #print self.r_vdw, self.gamma_w, self.alpha_ij, self.D_ij
        e_vdw = 0.0
        pow1 = math.pow((1/self.gamma_w), self.p_vdw1)
        pow2 = math.pow(r, self.p_vdw1)
        f13 = math.pow((pow1 + pow2), (1/self.p_vdw1))
        sum1 = 1 - f13/self.r_vdw
        exp1 = math.exp(self.alpha_ij * sum1)
        exp2 = math.exp(0.5 * self.alpha_ij * sum1)
        e_vdw = self.tap(r) * self.D_ij*(exp1 - 2*exp2)
        return e_vdw

    def coulomb(self, r):
        e_coulomb = 0.0
        a = math.pow(1/self.gamma_rij, 3)
        b = math.pow(r,3)
        c = math.pow(b + a, 1/3.0)
        tap = self.tap(r)
        e_coulomb = tap * self.q1 * self.q2 / c
        return e_coulomb

    def innerwall(self, r):
        e_inner = 0.0
        alpha = self.alpha_inner
        rin = self.r_inner
        D = self.D_inner
        if abs(rin) > 0.01:
            p1 = alpha * ( 1 - r/rin)
            e_inner = D * math.exp(p1)
        return e_inner

def bond_energy(bonds, ff):
    """assign the bond order parameters
    """
    for i in range(len(bonds)):
        b = bonds[i]
        # assign bond order parameters
        p = ff.eq2[i+1]
        b.a1 = p[0]
        b.a2 = p[1]
        b.r_o = float(p[2])
        b.p_bo1 = float(p[3])
        b.p_bo2 = float(p[4])
        b.r_pi = float(p[5])
        b.p_bo3 = float(p[6])
        b.p_bo4 = float(p[7])
        b.r_pipi = float(p[8])
        b.p_bo5 = float(p[9])
        b.p_bo6 = float(p[10])
        # assign bond energy parameters
        p = ff.eq6[i+1]
        b.D_e = float(p[2])
        b.p_be1 = float(p[3])
        b.p_be2 = float(p[4])
        b.D_ep = float(p[5])
        b.D_epp = float(p[6])
        # assign vdw parameters
        p = ff.eq23[i+1]
        b.rcut = 8.0
        b.D_ij = float(p[2])
        b.alpha_ij = float(p[3])
        b.r_vdw = float(p[4])
        b.r_vdw1 = 1.5591  # r_vdw1 is a globle parameter
        b.gamma_w = float(p[5])
        # assign charge parameters
        "to do"
        # assign inner wall
        p = ff.eq26[i+1]
        b.D_inner = float(p[2])
        b.alpha_inner = float(p[3])
        b.r_inner = float(p[4])

def plot_bo_be(r, bo, bo_s, bo_p, bo_pp, be, vdw, inner, name):
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(r, bo, '-o', label="BO")
    ax1.plot(r, bo_s, '-o', label="BO_sigma")
    if np.sum(bo_p) > 0:
        ax1.plot(r, bo_p, '-o', label="BO_pi")
    if np.sum(bo_pp) > 0:
        ax1.plot(r, bo_pp, '-o', label="BO_pipi")
    ax1.set_xlabel("Bond length ($\AA$)", size="x-large")
    ax1.set_ylabel("Bond order", size="x-large")
    ax1.set_title(name, size="x-large")
    ax1.legend()
    # plot bond order
    ax2.plot(r, be, '-o', label="be")
    #ax2.plot(r, vdw, '-o')
    ax2.plot(r, inner + vdw , '-o', label="inner")
    ax2.plot(r, vdw , '-o', label="vdw")
    #ax2.plot(r, be + vdw + inner, color="black", lw=3)
    ax2.plot(r, be + inner + vdw, color="black", lw=3)
    ax2.legend()
    # print min point
    tmp = be + inner + vdw
    print "re = %.3f"%r[tmp.argmin()],
    print "E = %.3f"%np.min(tmp)
    # save files
    for ii in range(100):
        fname = "bo.dat.%03d"%ii
        if not os.path.exists(fname):
            np.savetxt(fname, bo)
            break
    # plot bond order
    #plt.savefig("%s.eps"%name)
    #plt.savefig("%s.png"%name)
    plt.show()

def plot_bo_be2(r, bo, bo_s, bo_p, bo_pp, be, vdw, inner, name):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(r, bo, '-o', label="BO")
    ax1.set_xlabel("Bond length ($\AA$)", size="x-large")
    ax1.set_ylabel("Bond order", size="x-large")
    ax1.set_title(name, size="x-large")
    ax1.legend()
    # plot bond order
    ax2.plot(r, -be, '-o', label="be", color='red')
    ax2.set_ylim([min(-be), max(-be)])
    plt.show()

def plot_bond_order(bonds, r):
    """plot the bond order curve
    """

    counter = 1
    for i in range(len(bonds)):
        bo = []; bo_s = []; bo_p = []; bo_pp = []
        be = []; be_s = []; be_p = []; be_pp = []
        vdw = []
        inner = []
        if counter:
            for j in r:
                bonds[i].bond_order(j)
                bo.append(bonds[i].bo)
                bo_s.append(bonds[i].bo_s)
                bo_p.append(bonds[i].bo_p)
                bo_pp.append(bonds[i].bo_pp)
                bonds[i].bond_energy()
                be.append(bonds[i].be)
                be_s.append(bonds[i].be_s)
                be_p.append(bonds[i].be_p)
                be_pp.append(bonds[i].be_pp)
                vdw.append(bonds[i].vdw(j))
                inner.append(bonds[i].innerwall(j))
            bo = np.array(bo)
            bo_s = np.array(bo_s)
            bo_p = np.array(bo_p)
            bo_pp = np.array(bo_pp)
            be = np.array(be)
            be_s = np.array(be_s)
            be_p = np.array(be_p)
            be_pp = np.array(be_pp)
            vdw = np.array(vdw)
            inner = np.array(inner)
            name = bonds[i].a1 + "_" + bonds[i].a2
            plot_bo_be(r, bo, bo_s, bo_p, bo_pp, be, vdw, inner, name)
            #plot_bo_be2(r, bo, bo_s, bo_p, bo_pp, be, vdw, inner, name)
        counter += 1

def test():
    ff  = FfieldEq()
    bonds = []
    #init bond
    for i in range(len(ff.eq2) - 1):
        bond = Bond()
        bonds.append(bond)
    #assign bond order
    r = np.linspace(0.5, 2.5, 100)
    bond_energy(bonds, ff)
    plot_bond_order(bonds, r)


if __name__ == "__main__":
    test()

