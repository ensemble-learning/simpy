"""
call the eos program to calculate the bulk moldulus.
"""
import numpy as np
import math
import sys
import os

doc = """                         +---------------------------+
                         |     EOS Version 1.4.0     |
                         +---------------------------+

Equation of state (EOS) program for fitting energy-volume data. The following
variables are set in the file eos.in:

 cname               : name of crystal up to 256 characters
 natoms              : number of atoms in unit cell
 etype               : equation of state type (see below)
 vplt1, vplt2, nvplt : volume interval over which to plot energy, pressure etc.
                       as well as the number of points in the plot
 nevpt               : number of energy-volume points to be inputted
 vpt(i) ept(i)       : energy-volume points (atomic units)

Note that the input units are atomic - Bohr and Hartree (NOT Rydbergs).

The equations of state currently implemented are:
 1. Universal EOS (Vinet P et al., J. Phys.: Condens. Matter 1, p1941 (1989))
 2. Murnaghan EOS (Murnaghan F D, Am. J. Math. 49, p235 (1937))
 3. Birch-Murnaghan 3rd-order EOS (Birch F, Phys. Rev. 71, p809 (1947))
 4. Birch-Murnaghan 4th-order EOS
 5. Natural strain 3rd-order EOS (Poirier J-P and Tarantola A, Phys. Earth
    Planet Int. 109, p1 (1998))
 6. Natural strain 4th-order EOS
 7. Cubic polynomial in (V-V0)

--------------------------------------------------------------------------------
J. K. Dewhurst
August 2005
"""

A2B = math.pow(1.8897, 3) # angstrom to bohr
E2H = 0.037  # ev to hartree

x = np.loadtxt("vols")
y = np.loadtxt("results")

if not len(x) == len(y):
    sys.exit() 

x = x*A2B
y = y*E2H

o = open("eos.in", 'w')
o.write(' "eos calculation"     :cname\n')
o.write(" 1                            :natoms\n")
o.write(" 2                           : etype\n")
o.write(" %.2f  %.2f  %d  :vplt1, vplt2, nvplt\n"%(x[0], x[-1], len(x)))
o.write(" %d            :nevpt\n"%(len(y)))
for i in range(len(x)):
    o.write("%15.7f"%x[i])
    o.write("%15.7f"%y[i])
    o.write("\n")
o.close()

sig = os.system("eos")
if sig == 0:
    os.system("cat PARAM.OUT")
else:
    print "Error in calculation!"

