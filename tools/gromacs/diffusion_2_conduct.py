#@note: https://www.aqion.de/site/77
# ion       Di[m^2 s^-1]     A[S cm^2 mol^-1]
# H+        9.31e-9          349.6
# Na+       1.33e-9           50.0
# K+        1.96e-9           73.6
# OH-       5.27e-9          197.9
# Cl-       2.03e-9           76.2
# Br-       2.01e-9           75.5

#D: diffusion (10^-10 m^2/s)
#delta: d1, conductivity (Sm^2/mol or 10^4 S cm^2/mol)

L_2_m3 = 1e-3
F = 96485.33212 # Faraday constant, in C/mol
R = 8.314 # gas constant J/mol/K
NA = 6.02214076e23 #Avogadro constant

def test():
    D = 5.27e-9 # diffusion coefficients in 10^-10 m^2/s
    zi = 1.0 # ion valency
    F = 96485.33212 # Faraday constant, in C/mol
    R = 8.314 # gas constant J/mol/K
    T = 298.0 # in K
    RT = R*T

    #F^2/RT = 3.7554e6 s*S/mol

    #delta = (zi*ci*F^2*D)/(RT)

    d1 = zi*zi*F*F*D/RT # in Sm^2/mol

    print('molar limiting conductivity of ion i (S m^2/mol): %.4f'%d1)
    print('molar limiting conductivity of ion i (S cm^2/mol): %.2f'%(d1*1e4))

def diffusion_to_conductivity(D, zi, T, f, ci):
    RT = R*T
    d1 = zi*zi*F*F*D/RT # in Sm^2/mol
    d2 = d1*ci*f/L_2_m3 
    d2 = d2/100         # in S/cm

    print('ion concentration (M): %.2f'%(ci))
    print('dffusion coefficient of ion i (1e-5 cm^2/s): %.4f'%(D*1e9))
    print('molar limiting conductivity of ion i (S m^2/mol): %.4e'%d1)
    print('conductivity of ion i (S/cm): %.2e'%d2)
    print('conductivity of ion i (uS/cm): %.1f'%(d2*1e6))
    print('conductivity of ion i (mS/cm): %.2f'%(d2*1e3))

def single_ion():
    #@note: doi:10.1149/2.0461610jes
    """
    N = 1   # number of species in N
    L = 1.184 # the lenght of the cubic cell in nm
    D = 7.2e-11 # in m^2/s or 7.2e-7 cm^2/s
    zi = 1.0 # ion valency
    T = 298.0 # in K
    f = 0.64 # fraction of free ions
    ci = 0.52 # ion concentration, in mol/L
    diffusion_to_conductivity(D, zi, T, f, ci) 

    # DOI: 10.1149/2.0461610jes
    N = 300   # number of species in N
    L = 3.5608 # the lenght of the cubic cell in nm
    D = 0.072e-9 # diffusion coefficients in 10^-10 m^2/s

    zi = 1.0 # ion valency
    T = 298.0 # in K
    f = 0.64 # fraction of free ions
    """

    N = 300   # number of species in N
    L = 3.5608 # the lenght of the cubic cell in nm
    D = 0.0134e-9 # diffusion coefficients in 10^-10 m^2/s

    zi = 1.0 # ion valency
    T = 298.0 # in K
    f = 1.00 # fraction of free ions
    V = L*L*L # the volume of the system in nm^3
    ci = N/NA/V*1e24 # ion concentration, in mol/L

    if 0:
        ci = 0.52 # ion concentration, in mol/L

    diffusion_to_conductivity(D, zi, T, f, ci) 

if __name__ == "__main__":
    test()
    single_ion()
