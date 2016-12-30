"""
Plot fe.dat from VASP TI calculation
"""

import numpy as np
from scipy import integrate
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
    data = np.loadtxt("fe.dat") 
    data = data.transpose() 
    
    x, y, yerr = data[0], data[1], data[2]
    y_int = integrate.cumtrapz(y, x, initial=0)

    linestyle = mpltex.linestyle_generator(
                colors=["black", "green", "red"],
                lines=['-',':'], markers=['o','s'], 
                hollow_styles=[False, False, True, True],
                )

    fig, ax = plt.subplots(1)
    ax.plot(x, y_int, label="FE",
             **linestyle.next())

    ax.set_xlabel("r$_{C-O}$ ($\AA$)")
    ax.set_ylabel("FE (eV)")
    ax.legend(loc='best')

    ax2 = ax.twinx()
    ax2.errorbar(x, y, yerr=yerr, label="PMF", 
            **linestyle.next())
    ax2.set_ylabel("PMF (eV/$\AA$)", color="green")
    ax2.legend(loc='lower left')

    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)
    o = open("all.dat", "w")
    for i in range(len(x)):
        o.write("%12.4f%12.4f%12.4f%12.4f\n"%
                (x[i], y[i], yerr[i], y_int[i]))
    o.close()

plot_fig()


