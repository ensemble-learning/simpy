# -*- coding: utf-8 - *-
# https://pubs.acs.org/doi/abs/10.1021/acscatal.9b03239

import numpy as np
import os
from scipy.integrate import odeint

# Physical constants and conversion factors
# default is in eV
J2eV = 6.24150974E18 # eV/J
Na = 6.0221415E23 # mol-1
h = 6.626068E-34 * J2eV # in eV*s
kb = 1.3806503E-23 * J2eV # in eV/K

# Reaction conditions
T = 300.0
kbT = kb * T # in eV

class KMC():
    def __init__(self,):
        self.species = []
        self.is_surface = []
        self.species_sur = []
        self.c = {}

        self.rxns = []
        self.dErxn = []
        self.dSrxn = []
        self.dGrxn = []
        self.dEact = []
        self.dSact = []
        self.dGact = []
        self.kf = []
        self.kb = []

def read_rxn(kmc):
    f = open('rxn.csv', 'r')
    for i in f:
        if not i.strip().startswith('#'):
            tokens = i.strip().split(',')
            kmc.dErxn.append(float(tokens[1].strip()))
            kmc.dSrxn.append(float(tokens[2].strip()))
            kmc.dEact.append(float(tokens[3].strip()))
            kmc.dSact.append(float(tokens[4].strip()))
            r = tokens[0].split('<-->')[0]
            p = tokens[0].split('<-->')[1]
            reactants = [ii.strip() for ii in r.split('+')]
            for ii in reactants:
                if not ii == '*' and not ii in kmc.species:
                        kmc.species.append(ii)
            products = [ii.strip() for ii in p.split('+')]
            for ii in products:
                if not ii == '*'and not ii in kmc.species:
                        kmc.species.append(ii)
            kmc.rxns.append([reactants, products])
    f.close()

    for n, i in enumerate(kmc.species):
        tokens = i.strip()
        if '*' in tokens:
            kmc.species_sur.append(tokens)
            kmc.is_surface.append(1)
        else:
            kmc.is_surface.append(0)
        kmc.c[tokens] = 0.0

    for e, s in zip(kmc.dErxn, kmc.dSrxn):
        kmc.dGrxn.append(e-T*s)
    for e, s in zip(kmc.dEact, kmc.dSact):
        kmc.dGact.append(e-T*s)
    # write species.dat
    o = open('species.template', 'w')
    for i in kmc.species:
        o.write('%10s %10.4f\n'%(i, 0.))
    o.close()

def get_rate_constants(kmc):
    # Calculate equilibrium and rate constants
    K = []
    for i in range(len(kmc.dGrxn)):
        K.append(np.exp(-kmc.dGrxn [i]/kbT))
        kmc.kf.append(kbT/h * np.exp(-kmc.dGact [i]/kbT))
    for i in range(len(kmc.dGrxn)):
        kmc.kb.append(kmc.kf[i]/K[i]) # enforce thermodynamic consistency

def assign_initial_c(kmc):
    if os.path.exists('species.dat'):
        f = open('species.dat', 'r')
        for i in f:
            tokens = i.strip().split()
            kmc.c[tokens[0]] = float(tokens[1])

def get_rates(theta,kmc):
    # returns the rates depending on the current coverages theta

    tstar = 1.0
    for n, i in enumerate(kmc.species_sur):
        kmc.c[i] = theta[n]
        tstar = tstar - kmc.c[i]
    kmc.c['*'] = tstar

    # Caluclate the rates:
    rate = [0.]*len(kmc.rxns) # array with 21 zeros, one for each reaction
   
    for n, rxn in enumerate(kmc.rxns):
        rf, rb = kmc.kf[n], kmc.kb[n]
        for r in rxn[0]:
            rf = rf*kmc.c[r]
        for p in rxn[1]:
            rb = rb*kmc.c[p]
        rate[n] = rf - rb
    return rate

def get_odes(theta,t,kmc):
    # returns the system of ODEs d(theta)/dt
    # calculated at the current value of theta (and time t, not used)
    rate = get_rates(theta,kmc) # calculate the current rates
    # Time derivatives of theta
    dt = [0]*14 # 14 surface species
   
    for ni, i in enumerate(kmc.species_sur):
        for nj, j in enumerate(kmc.rxns): 
            for ii in j[0]:
                if ii == i:
                    dt[ni] += -rate[nj]
            for ii in j[1]:
                if ii == i:
                    dt[ni] += rate[nj]
    return dt

def print_output(theta, kmc):
    # Prints the solution of the model
    get_rate_constants(kmc)
    rates = get_rates(theta,kmc)
   
    for n,rate in enumerate(rates):
        print('%-20s <--> %20s'%(' + '.join(kmc.rxns[n][0]), ' + '.join(kmc.rxns[n][1])),
              ": rate=%.4e, kf=%.4e, kb=%.4e"%(rate,kmc.kf[n], kmc.kb[n]))

    print("The coverages are:")
    for n, i in enumerate(kmc.species_sur):
        print('%10s'%i, ': ', '%.4f'%theta[n])

if __name__ == '__main__':
    kmc = KMC()

    read_rxn(kmc)
    assign_initial_c(kmc)
    get_rate_constants(kmc) 

    theta0 = [0.]*len(kmc.species_sur)
    # Integrate the ODEs for 1E5 sec (enough to reach steady-state)
    theta = odeint(get_odes, # system of ODEs 
                   theta0, # initial guess 
                   [0,1E5 ], # Integrate ODE until we reach time of interest (t') 
                   args = (kmc,), # additional arguments to get_odes() 
                   h0 = 1E-36, # initial time step 
                   mxstep =10000, # maximum number of steps
                   rtol = 1E-12, # relative tolerance
                   atol = 1E-15) # absolute tolerance
    theta_final = theta[-1,:]
    print_output(theta_final, kmc)

