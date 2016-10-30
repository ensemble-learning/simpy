"""
Plot radial distribution function in acs style.
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
    fig, ax = plt.subplots(1)
    linestyle = mpltex.linestyle_generator(
                colors=["black", "black", "red", "red"],
                lines=['-',':'], markers=['o','o','^', '^'], 
                hollow_styles=[False, False, True, True],
                )

    # plot data 1
    data = np.loadtxt("gofr_1.dat")
    data = data.transpose()
    x, y = [], []
    for i in range(len(data[0])):
        if i <= 50:
            x.append(data[0][i])
            y.append(data[1][i])
        else:
            if i % 10 == 0:
                x.append(data[0][i])
                y.append(data[1][i])
        
    ax.plot(x, y, label="I$^-$(bulk)-H",
            **linestyle.next())

    ax2 = ax.twinx()
    ax2.plot(data[0][::10], data[2][::10],
             **linestyle.next())

    # plot data 2
    data = np.loadtxt("gofr_2.dat")
    data = data.transpose()
    x, y = [], []
    for i in range(len(data[0])):
        if i <= 50:
            x.append(data[0][i])
            y.append(data[1][i])
        else:
            if i % 10 == 0:
                x.append(data[0][i])
                y.append(data[1][i])
    ax.plot(x, y, label="I*(interface)-H", 
            **linestyle.next())

    # add annodate
    ax2.plot(data[0][::5], data[2][::5],
             **linestyle.next())

    ax2.plot([3.0,8.0], [5.0, 5.0], ls="dotted",
             lw=0.5, color="black")

    ax2.text(7.0, 5.2, "N = 5")

    # set up axis
    ax.set_xlabel("r$_{I-H}$ ($\AA$)")
    ax.yaxis.set_ticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_ylabel("g(r)")
    ax.legend(loc='upper right', bbox_to_anchor=(0.85, 0.95))

    ax2.set_ylim([0.0, 8.0])
    ax2.set_xlim([1.0, 8.0])
    ax2.set_ylabel("$\int$g(r) (N)") 

    # output
    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)

plot_fig()
