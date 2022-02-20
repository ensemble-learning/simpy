# -*- coding: utf-8 - *-
# https://pubs.acs.org/doi/abs/10.1021/acscatal.9b03239

import numpy as np
import os

class KMC():
    def __init__(self,):
        self.species = []
        self.is_surface = []
        self.species_sur = []
        self.dErxn = []
        self.dSrxn = []
        self.dGrxn = []
        self.dEact = []
        self.dSact = []
        self.dGact = []
        self.rxns = []
        self.c = {}

#To run the model at any pH, remove the comment at the start of the line in both reaction energies and activation barriers sections
#Current settings are at pH = 4 and T = 300 K

# Reaction conditions
T = 300.0
A = 0.8 #Scaling parameter for activation barriers

# Physical constants and conversion factors
J2eV = 6.24150974E18 # eV/J
Na = 6.0221415E23 # mol-1
h = 6.626068E-34 * J2eV # in eV*s
kb = 1.3806503E-23 * J2eV # in eV/K
kbT = kb * T # in eV

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
    if 0:
        for i in kmc.species:
            print(i, '0.00')

def get_rate_constants(kmc):
    # Calculate equilibrium and rate constants
    K, kf, kr = [], [], []
    for i in range(len(kmc.dGrxn)):
        K.append(np.exp(-kmc.dGrxn [i]/kbT))
        kf.append(kbT/h * np.exp(-kmc.dGact [i]/kbT))
    for i in range(len(kmc.dGrxn)):
        kr.append(kf[i]/K[i]) # enforce thermodynamic consistency
    return kf, kr

def assign_initial_c(kmc):
    if os.path.exists('species.dat'):
        f = open('species.dat', 'r')
        for i in f:
            tokens = i.strip().split()
            kmc.c[tokens[0]] = float(tokens[1])

kmc = KMC()
read_rxn(kmc)
assign_initial_c(kmc)

def get_rates(theta,kmc):
    # returns the rates depending on the current coverages theta

    nc = 0
    tstar = 1.0
    for n, i in enumerate(kmc.species_sur):
        kmc.c[i] = theta[nc]
        tstar = tstar - kmc.c[i]
        nc += 1
    kmc.c['*'] = tstar

    # Caluclate the rates:
    rate = [0]*len(kmc.rxns) # array with 21 zeros, one for each reaction
   
    for n, rxn in enumerate(kmc.rxns):
        rf, rb = 1.0, 1.0
        for r in rxn[0]:
            rf = rf*kmc.c[r]
        for p in rxn[1]:
            rb = rb*kmc.c[p]
        rate[n] = kf[n]*rf - kr[n]*rb
        if n == 2:
            print(rf, rb, rxn)
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
                    dt[ni] = dt[ni] - rate[nj]
            for ii in j[1]:
                if ii == i:
                    dt[ni] = dt[ni] + rate[nj]
    return dt

def print_output(theta, kmc, kf, kr):
    # Prints the solution of the model
    get_rate_constants(kmc)
    rates = get_rates(theta,kmc)
   
    print
    for r,rate in enumerate(rates):
        print ("Step",r,": rate =",rate,", kf =",kf[r],", kr=",kr[r])
    print
    print ("The coverages for NO*, OH*, O*, N*, NH*, NH2*, NH3*, H2O*, H*, NO2*, NOH*, HNO*, HNOH*, N2H4* are:")
    for t in theta:
        print (t)
    print
    print ("Rate of reactions")
    for r, rate in enumerate(rates):
   	    print (rate)

#Reaction energies
#To run the model at any pH, remove the comment at the start of the line in both reaction energies and activation barriers sections

(kf, kr) = get_rate_constants(kmc) 

# Use scipy’s odeint to solve the system of ODEs
from scipy.integrate import odeint
# As initial guess we assume an empty surface
theta0 = (0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0.)
# Integrate the ODEs for 1E5 sec (enough to reach steady-state)
theta = odeint(get_odes, # system of ODEs 
               theta0, # initial guess 
               [0,1E5 ], # Integrate ODE until we reach time of interest (t') 
               args = (kmc,), # additional arguments to get_odes() 
               h0 = 1E-36, # initial time step 
               mxstep =1000, # maximum number of steps
               rtol = 1E-12, # relative tolerance
               atol = 1E-15) # absolute tolerance
print_output(theta [-1,:], kmc, kf, kr)
theta_final = theta [-1,:]

