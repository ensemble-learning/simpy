import math
import numpy as np
from scipy import interpolate

def getY(start, end, n):
    interp = []

    a = np.linspace(0.8, 1.2, 11)
    x = a*a*a
    data = np.loadtxt("results")
    data = data.transpose()
    y = data

    tck = interpolate.UnivariateSpline(x, y, w=None, bbox=[None, None], k=5, s=2)

    for i in np.linspace(start, end, n):
        interp.append(tck(i))

    return interp

def getV():
    f = open("./scan_05/CONTCAR", 'r')
    f.readline()
    scale = f.readline()
    a = f.readline().strip().split()
    b = f.readline().strip().split()
    c = f.readline().strip().split()
    v = float(a[0]) * float(b[1]) * float(c[2])
    return v


A2B = math.pow(1.8897, 3)
E2H = 0.037

V0 = getV()

start = 0.9
end = 1.1

npoints = 1000

n = 20

x = np.linspace(start, end, n)
x = x*V0*A2B
y = getY(start, end, n)
y = np.array(y)
y = y*E2H

o = open("eos.in", 'w')
o.write(' "genenerated from python"     :cname\n')
o.write(" 28                            :natoms\n")
o.write(" 2                           : etype\n")
o.write(" %.2f  %.2f     %d    :vplt1, vplt2, nvplt\n"%(x[0], x[-1], npoints))
o.write(" %d            :nevpt\n"%(len(y)))

INP = """ "Silicon"                    : cname
 28                            : natoms
 2                            : etype
 1 2 1000             : vplt1, vplt2, nvplt
 11                            : nevpt
"""

for i in range(n):
    o.write("%15.7f"%x[i])
    o.write("%15.7f"%y[i])
    o.write("\n")
o.close()

