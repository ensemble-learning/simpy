#!/usr/bin/env python
#
def to_plot(hw,ab,gam=0.001,type='Gaussian'):
    import numpy as np
    #
    fmin = min(hw)
    fmax = max(hw)
    erange = np.arange(fmin-40*gam,fmax+40*gam,gam/10)#np.arange(fmin-40*gam,fmax+40*gam,gam/10)
    spectrum = 0.0*erange
    for i in range(len(hw)):
        if type=='Gaussian':
            spectrum += (2*np.pi)**(-.5)/gam*np.exp(np.clip(-1.0*(hw[i]-erange)**2/(2*gam**2),-300,300))
        elif type=='Lorentzian':
            spectrum += ab[i]*1/np.pi*gam/((hw[i]-erange)**2+gam**2)
    #
    return erange, spectrum
#
if __name__ == '__main__':
    import numpy as np
    import sys
    #
    hw=np.genfromtxt(sys.argv[1], dtype=float)
    cm1 = hw[:,1]
    int1 = hw[:,4]
    int1 /= np.max(np.abs(int1),axis=0)
    Es1,Spectrum1 = to_plot(cm1, int1, 15.0, 'Lorentzian')
#
    filename = sys.argv[1]+"-broaden.dat"
    print "Writing", filename
    f = open(filename,'w')
    f.write('# freq/cm-1  Intensity \n')
    for i in range(len(Es1)):
        f.write('%.5e   %.5e\n' % (Es1[i],Spectrum1[i]))
    f.close()
