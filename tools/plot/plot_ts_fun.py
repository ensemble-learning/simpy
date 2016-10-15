import sys
import numpy as np
#import matplotlib.pyplot as plt

def plot_rxn(x0, x1, e_t0, e_ts, e_t1):

    n_points = 30
    a = 0.5*(x0 + x1)
    b = e_t0 + e_ts

    k1 = (b-e_t0)/((x0-a)*(x0-a))
    print "fun 1: y = -%.3f*(x-%.3f)^2 + %.3f"%(k1, a, b)

    k2 = (b-e_t1)/((x1-a)*(x1-a))
    print "fun 2: y = -%.3f*(x-%.3f)^2 + %.3f"%(k2, a, b)

    #x = np.linspace(x0,a,n_points)
    #y = -k1*(x - a) * (x - a) + b          
    #plt.plot(x,y, color="black") 

    #x = np.linspace(a,x1,n_points)
    #y = -k2*(x - a) * (x - a) + b          
    #plt.plot(x,y, color="black") 

    #plt.scatter([x0, a, x1], [e_t0, e_ts, e_t1])

    #plt.show()

print "input x0, x1, e_t0, e_ts, e_t1"
if len(sys.argv) == 6:
    x0 = float(sys.argv[1])
    x1 = float(sys.argv[2])
    e_t0 = float(sys.argv[3])
    e_ts = float(sys.argv[4])
    e_t1 = float(sys.argv[5])

    plot_rxn(x0, x1, e_t0, e_ts, e_t1)

