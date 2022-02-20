# -*- coding: utf-8 - *-
# https://pubs.acs.org/doi/abs/10.1021/acscatal.9b03239

import numpy as np
import os

class KMC():
    def __init__(self,):
        self.species = []
        self.is_surface = []
        self.c = {} # concentrations
        self.tstar = 0.0
        self.dErxn = []
        self.dSrxn = []
        self.dGrxn = []
        self.dEact = []
        self.dSact = []
        self.dGact = []
        self.rxns = []
        self.rates = []
        self.kf = []
        self.kr = []
        self.K = []

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

PHNO2 = 0.03 # bar
PH2 = 0.02 # bar
PH2O = 1 # bar
PNH3 = 0.0 #bar
PN2H4 = 0.0 #bar

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
    for i in range(len(kmc.dGrxn)):
        kmc.K.append(np.exp(-kmc.dGrxn [i]/kbT))
        kmc.kf.append(kbT/h * np.exp(-kmc.dGact [i]/kbT))
    for i in range(len(kmc.dGrxn)):
        kmc.kr.append(kmc.kf[i]/kmc.K[i]) # enforce thermodynamic consistency

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
    # Extract elements of theta and assign them
    # to more meaningful variables

    nc = 0
    tstar = 1.0
    for n, i in enumerate(kmc.species):
        if kmc.is_surface[n] == 1:
            kmc.c[i] = theta[nc]
            tstar = tstar - kmc.c[i]
            nc += 1
    kmc.c['*'] = tstar

    tNO = theta [0] # theta of NO
    tOH = theta [1] # theta of OH
    tO = theta [2] # theta of O
    tN = theta [3] # theta of N
    tNH = theta [4] # theta of NH
    tNH2 = theta [5] # theta of NH2
    tNH3 = theta [6] # theta of NH3
    tH2O = theta [7] # theta of H2O
    tH = theta [8] # theta of H
    tNO2 = theta [9] # theta of NO2
    tNOH = theta [10] #theta of NOH
    tHNO = theta [11] #theta of HNO
    tHNOH = theta [12] #theta of HNOH
    tN2H4 = theta [13] #theta of N2H4
    tstar = 1.0 - tNO - tOH - tO - tN - tNH - tNH2 - tNH3 - tH2O - tH - tNO2 - tNOH - tHNO - tHNOH - tN2H4 # site balance for tstar
        
    # Caluclate the rates:
    rate = [0]*len(kmc.rxns) # array with 21 zeros, one for each reaction
   
    for n, rxn in enumerate(kmc.rxns):
        rf, rb = 1.0, 1.0
        for r in rxn[0]:
            rf = rf*kmc.c[r]
        for p in rxn[1]:
            rb = rb*kmc.c[p]
        rate[n] = kmc.kf[n]*rf - kmc.kr[n]*rb
    #print('before: ', rate)

    kf = kmc.kf
    kr = kmc.kr
    # Caluclate the rates:
    rate = [0]*21 # array with 21 zeros, one for each reaction

    rate [0] = kf [0] * PHNO2 * tstar * tstar   - kr [0] * tNO * tOH
    rate [1] = kf [1] * PH2 * tstar * tstar     - kr [1] * tH * tH
    rate [2] = kf [2] * tOH * tH                - kr [2] * tH2O * tstar
    rate [3] = kf [3] * tNO * tstar             - kr [3] * tN * tO
    rate [4] = kf [4] * tO * tH                 - kr [4] * tOH * tstar
    rate [5] = kf [5] * tN * tH                 - kr [5] * tNH * tstar
    rate [6] = kf [6] * tNH * tH                - kr [6] * tNH2 * tstar
    rate [7] = kf [7] * tNH2 * tH               - kr [7] * tNH3 * tstar
    rate [8] = kf [8] * tNH3                    - kr [8] * PNH3 * tstar
    rate [9] = kf [9] * tH2O                    - kr [9] * PH2O * tstar
    rate [10] = kf [10] * PHNO2 * tstar * tstar - kr [10] * tNO2 * tH
    rate [11] = kf [11] * tNO2 * tstar          - kr [11] * tNO * tO
    rate [12] = kf [12] * tNO * tH              - kr [12] * tNOH * tstar
    rate [13] = kf [13] * tNO * tH              - kr [13] * tHNO * tstar
    rate [14] = kf [14] * tNOH * tH             - kr [14] * tHNOH * tstar
    rate [15] = kf [15] * tNOH * tstar          - kr [15] * tN * tOH
    rate [16] = kf [16] * tHNO * tH             - kr [16] * tHNOH * tstar
    rate [17] = kf [17] * tHNO * tstar          - kr [17] * tNH * tO
    rate [18] = kf [18] * tHNOH * tstar         - kr [18] * tNH * tOH
    rate [19] = kf [19] * tNH2 * tNH2           - kr [19] * tN2H4 * tstar
    rate [20] = kf [20] * tN2H4                 - kr [20] * PN2H4 * tstar
    #print('aafter: ',rate)

    return rate

def get_odes(theta,t,kmc):
    # returns the system of ODEs d(theta)/dt
    # calculated at the current value of theta (and time t, not used)
    rate = get_rates(theta,kmc) # calculate the current rates
    # Time derivatives of theta
    dt = [0]*14 # 14 surface species
   
    dt [0] = rate [0] - rate [3] + rate [11] - rate [12] - rate [13]  # d(tNO)/dt
    dt [1] = rate [0] + rate [15] + rate [18] - rate [2] + rate [4]   # d(tOH)/dt
    dt [2] = rate [3] + rate [17] + rate [11] - rate [4]              # d(tO)/dt
    dt [3] = rate [3] + rate [15] - rate [5]                          # d(tN)/dt
    dt [4] = rate [5] + rate [17] + rate [18] - rate [6]              # d(tNH)/dt
    dt [5] = rate [6] - rate [7] - 2 * rate [19]                      # d(tNH2)/dt
    dt [6] = rate [7] - rate [8]                                      # d(tNH3)/dt
    dt [7] = rate [2] - rate [9]                                      # d(tH2O)/dt
    dt [8] = 2 * rate [1] + rate [10] - rate [2] - rate [4] - rate [5] - rate [6] - rate [7] - rate [12] - rate[13] - rate [14] - rate [16]  # d(tH)/dt
    dt [9] = rate [10] - rate [11]                                    # d(tNO2)/dt
    dt [10] = rate [12] - rate [14] - rate [15]                       # d(NOH)/dt 
    dt [11] = rate [13] - rate [16] - rate [17]                       # d(HNO)/dt
    dt [12] = rate [14] + rate [16] - rate [18]                       # d(HNOH)/dt
    dt [13] = rate [19] - rate [20]                                   # d(N2H4)/dt 
    return dt

def print_output(theta, kmc):
    # Prints the solution of the model
    get_rate_constants(kmc)
    rates = get_rates(theta,kmc)
   
    print
    for r,rate in enumerate(rates):
        print ("Step",r,": rate =",rate,", kf =",kmc.kf[r],", kr=",kmc.kr[r])
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

def free_energies(kmc):
    E = [-2.864, -1.165, 0.06, -0.85, 0.22, -0.256, 0.323, -0.153, 0.572, 0.656, -1.734, -1.64, 0.501, 0.599, 0.449, -1.133, 0.351, -1.708, -1.839, 1.039, 1.449] #4
    #calculate the Gibbs free energy for each reaction using Erxn and Srxn
    S = [-1.965E-3, -1.836E-3, 2.737E-4, -1.336E-4, 9.791E-5, -4.539E-7, 7.032E-5, 2.258E-4, 1.533E-3, 1.4021E-3, -2.954E-3, -8.705E-5, 1.268E-4, -7.899E-5, -1.26E-4, -1.625E-4, 7.978E-5, -5.505E-5, -3.693E-5, 6.272E-5, 1.932E-3]
    kmc.Grxn = [ E [i] - T * S [i] for i in range (len(E))]

free_energies(kmc)


def Activation_free_energies():
    #Activation Barriers
    E = [0.045, 0.000, 0.825 *A , 1.432 *A , 1.250 *A , 0.892 *A , 1.164 *A , 1.053 *A , 0.572, 0.656, 0.0, 0.267*A, 1.621 *A, 1.492*A, 1.34*A, 0.70*A, 0.73*A, 0.75*A, 0.0, 1.524 *A, 1.449 ] #pH=4
    S = [-5.133E-4, -3.776E-4, 8.264E-5, -1.107E-4, 8.137E-5, 4.401E-5, 8.583E-5, 8.480E-5, 1.063E-3, 9.301E-4, -1.397E-3, -6.158E-5, -2.03E-5, 3.36E-5 , 6.285E-5, -1.68E-4, -2.29E-4, -5E-5, 1E-8, 4.47E-5 , 1.349E-3]
    Gact = [ E [i]  - T * S [i]  for i in range (len(E))]
    return Gact
   
dGact = Activation_free_energies()
print (dGact)
        
get_rate_constants(kmc) 

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
print_output(theta [-1,:], kmc)
theta_final = theta [-1,:]

