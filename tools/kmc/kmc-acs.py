# -*- coding: utf-8 - *-
# https://pubs.acs.org/doi/abs/10.1021/acscatal.9b03239

import numpy as np

#To run the model at any pH, remove the comment at the start of the line in both reaction energies and activation barriers sections
#Current settings are at pH = 4 and T = 300 K

# Reaction conditions
T = 300 # K
A = 1.0  #Scaling parameter for activation barriers

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

def get_rate_constants(dGrxn,dGact):
    # Calculate equilibrium and rate constants
    K =[0]*len(dGrxn) # equilibrium constants for 21 reactions
    kf = [0]*len(dGact) # forward rate constants for 21 reactions
    kr = [0]*len(dGact) # reverse rate constants for 21 reactions
   
    for i in range(len(dGrxn)):
        K [i] = np.exp(-dGrxn [i]/kbT)
        kf [i] = kbT/h * np.exp(-dGact [i]/kbT)
        kr [i] = kf [i]/K [i] # enforce thermodynamic consistency
    return (kf,kr)

def get_rates(theta,kf,kr):
    # returns the rates depending on the current coverages theta
    # Extract elements of theta and assign them
    # to more meaningful variables
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
    rate = [0]*21 # array with 21 zeros, one for each reaction
   
    # HNO2 + 2* <--> NO* + OH*
    rate [0] = kf [0] * PHNO2 * tstar * tstar   - kr [0] * tNO * tOH        
    # H2 + 2* <--> 2H*
    rate [1] = kf [1] * PH2 * tstar * tstar     - kr [1] * tH * tH          
    # OH* + H* <--> H2O* + * 
    rate [2] = kf [2] * tOH * tH                - kr [2] * tH2O * tstar     
    # NO* + * <--> N* + O*
    rate [3] = kf [3] * tNO * tstar             - kr [3] * tN * tO          
    # O* + H* <--> OH* + *
    rate [4] = kf [4] * tO * tH                 - kr [4] * tOH * tstar      
    # N* + H* <--> NH* + *
    rate [5] = kf [5] * tN * tH                 - kr [5] * tNH * tstar      
    # NH* + H* <--> NH2* + *
    rate [6] = kf [6] * tNH * tH                - kr [6] * tNH2 * tstar     
    # NH2* + H* <--> NH3* + * 
    rate [7] = kf [7] * tNH2 * tH               - kr [7] * tNH3 * tstar     
    # NH3* <--> NH3 + *
    rate [8] = kf [8] * tNH3                    - kr [8] * PNH3 * tstar     
    # H2O* <--> H2O + *
    rate [9] = kf [9] * tH2O                    - kr [9] * PH2O * tstar     
    # HNO2 + 2* <--> NO2* + H*
    rate [10] = kf [10] * PHNO2 * tstar * tstar - kr [10] * tNO2 * tH       
    # NO2* + * <--> NO* + O*
    rate [11] = kf [11] * tNO2 * tstar          - kr [11] * tNO * tO        
    # NO* + H* <--> NOH*
    rate [12] = kf [12] * tNO * tH              - kr [12] * tNOH * tstar    
    # NO* + H* <--> HNO*
    rate [13] = kf [13] * tNO * tH              - kr [13] * tHNO * tstar    
    # NOH* + H* <--> HNOH*
    rate [14] = kf [14] * tNOH * tH             - kr [14] * tHNOH * tstar   
    # NOH* + * <--> N* + OH*
    rate [15] = kf [15] * tNOH * tstar          - kr [15] * tN * tOH        
    # HNO* + H* <--> HNOH*
    rate [16] = kf [16] * tHNO * tH             - kr [16] * tHNOH * tstar   
    # HNO* + * <--> NH* + O*
    rate [17] = kf [17] * tHNO * tstar          - kr [17] * tNH * tO        
    # HNOH* + * <--> NH* + OH*
    rate [18] = kf [18] * tHNOH * tstar         - kr [18] * tNH * tOH       
    # 2NH2* <--> N2H4* + *
    rate [19] = kf [19] * tNH2 * tNH2           - kr [19] * tN2H4 * tstar   
    # N2H4* <--> N2H4 + *
    rate [20] = kf [20] * tN2H4                 - kr [20] * PN2H4 * tstar   

    return rate

def get_odes(theta,t,kf,kr):
    # returns the system of ODEs d(theta)/dt
    # calculated at the current value of theta (and time t, not used)
    rate = get_rates(theta,kf,kr) # calculate the current rates
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

def print_output(theta):
    # Prints the solution of the model
    (kf,kr) = get_rate_constants(dGrxn,dGact)
    rates = get_rates(theta,kf,kr)
   
    print
    for r,rate in enumerate(rates):
        print ("Step",r,": rate =",rate,", kf =",kf [r],", kr=",kr [r])
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

def free_energies():
    E = [-2.864, -1.165, 0.06, -0.85, 0.22, -0.256, 0.323, -0.153, 0.572, 0.656, -1.734, -1.64, 0.501, 0.599, 0.449, -1.133, 0.351, -1.708, -1.839, 1.039, 1.449] #4
    #calculate the Gibbs free energy for each reaction using Erxn and Srxn
    S = [-1.965E-3, -1.836E-3, 2.737E-4, -1.336E-4, 9.791E-5, -4.539E-7, 7.032E-5, 2.258E-4, 1.533E-3, 1.4021E-3, -2.954E-3, -8.705E-5, 1.268E-4, -7.899E-5, -1.26E-4, -1.625E-4, 7.978E-5, -5.505E-5, -3.693E-5, 6.272E-5, 1.932E-3]
    Grxn = [ E [i] - T * S [i] for i in range (len(E))]
    return Grxn

dGrxn = free_energies()
print (dGrxn)


def Activation_free_energies():
    #Activation Barriers
    E = [0.045, 0.000, 0.825 *A , 1.432 *A , 1.250 *A , 0.892 *A , 1.164 *A , 1.053 *A , 0.572, 0.656, 0.0, 0.267*A, 1.621 *A, 1.492*A, 1.34*A, 0.70*A, 0.73*A, 0.75*A, 0.0, 1.524 *A, 1.449 ] #pH=4
    S = [-5.133E-4, -3.776E-4, 8.264E-5, -1.107E-4, 8.137E-5, 4.401E-5, 8.583E-5, 8.480E-5, 1.063E-3, 9.301E-4, -1.397E-3, -6.158E-5, -2.03E-5, 3.36E-5 , 6.285E-5, -1.68E-4, -2.29E-4, -5E-5, 1E-8, 4.47E-5 , 1.349E-3]
    Gact = [ E [i]  - T * S [i]  for i in range (len(E))]
    return Gact
   
dGact = Activation_free_energies()
print (dGact)
        
(kf,kr) = get_rate_constants(dGrxn,dGact) 

# Use scipy’s odeint to solve the system of ODEs
from scipy.integrate import odeint
# As initial guess we assume an empty surface
theta0 = (0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0.)
# Integrate the ODEs for 1E5 sec (enough to reach steady-state)
theta = odeint(get_odes, # system of ODEs 
               theta0, # initial guess 
               [0,1E5 ], # Integrate ODE until we reach time of interest (t') 
               args = (kf,kr), # additional arguments to get_odes() 
               h0 = 1E-36, # initial time step 
               mxstep =10000, # maximum number of steps
               rtol = 1E-12, # relative tolerance
               atol = 1E-15) # absolute tolerance
print_output(theta [-1,:])
theta_final = theta [-1,:]

'''
# Degree of rate control
print("Degree of Rate Control")

#Degree of Rate Control, Xrc - run the mkm model for a short timestep picking up where the main mkm stopped
theta_o = odeint(get_odes, # system of ODEs
                 theta_final, # state of the surface at t'
                 [0,1], # short time span
                 args = (kf,kr), # additional arguments to get_odes()
                 h0 = 1E-36, # initial time step
                 mxstep =9000000, # maximum number of steps
                 rtol = 1E-12, # relative tolerance
                 atol = 1E-15) # absolute tolerance

ro = get_rates(theta[-1,:],kf,kr)[8] #save rate
ratesdrc = [0]*21 #initialize DRC rates for each Step
Xrc = [0]*21 #initialize Xrc for each Step
x = 0.002 #set change in barrier height


for s in range(len(dGact)):
    dGactdrc = dGact[:] #reset barriers
    dGactdrc[s]=dGact[s]-x #modify barrier of step "s"
    (kfdrc,krdrc) = get_rate_constants(dGrxn,dGactdrc) #get rate constants with modified barrier for step "s"


    thetadrc = odeint(get_odes, # system of ODEs
                      theta_final, # state of the surface at t'
                      [0,1], # short time span
                      args = (kfdrc,krdrc), # additional arguments to get_odes()
                      h0 = 1E-36, # initial time step
                      mxstep =9000000, # maximum number of steps
                      rtol = 1E-12, # relative tolerance
                      atol = 1E-15) # absolute tolerance

    ratesdrc[s] = get_rates(thetadrc[-1,:],kfdrc,krdrc)[8] #compute new rate 
    Xrc[s] = ((ratesdrc[s] - ro) / ro)*(kf[s] / (kfdrc[s] - kf[s])) #compute Xrc for step "s"
    print ("Step",s,": Xrc =",Xrc[s])

print ("Sum:", sum(Xrc))

# Degree of transient rate control
print("Degree of Transient Rate Control")
ro = get_rates(theta_final,kf,kr)[8] #save rate
ratesdrc = [0]*21 #initialize DRC rates for each Step
Xrc = [0]*21 #initialize Xrc for each Step
x = 0.002 #set change in barrier height


for s in range(len(dGact)):
    dGactdrc = dGact[:] #reset barriers
    dGactdrc[s]=dGact[s]-x #modify barrier of Step "s"
    (kfdrc,krdrc) = get_rate_constants(dGrxn,dGactdrc) #get rate constants with modified barrier for step "s"


    thetadrc = odeint(get_odes, # system of ODEs
                      theta0, # state of the surface at t = 0
                      [0,1E5], # time span t = 0 to t'
                      args = (kfdrc,krdrc), # additional arguments to get_odes()
                      h0 = 1E-36, # initial time step
                      mxstep =9000000, # maximum number of steps
                      rtol = 1E-12, # relative tolerance
                      atol = 1E-15) # absolute tolerance

    ratesdrc[s] = get_rates(thetadrc[-1,:],kfdrc,krdrc)[8] #compute new rate 
    Xrc[s] = ((ratesdrc[s] - ro) / ro)*(kf[s] / (kfdrc[s] - kf[s])) #compute Xrc for step "s"
    print ("Step",s,": Xrc, transient =",Xrc[s])

print ("Sum:", sum(Xrc))
'''
