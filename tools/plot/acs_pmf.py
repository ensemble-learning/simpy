import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
    data = np.loadtxt("all.dat", 
                    dtype={'names': ('folder', 'x', 'dA/dx', 'err', 'A'),
                         'formats': ('S20', 'f4', 'f4', 'f4', 'f4')})
    xcut = 27.5 # free energy beyong 27.5 is too large
    x0 = 28.36  # the top layer of Cu
    
    x, y = [], []
    for i in range(len(data)):
        tx = data[i][1]
        ty = data[i][4]
        if tx < xcut:
            x.append(x0-tx)
            y.append(data[i][4])

    fig, ax = plt.subplots(1)
    linestyle = mpltex.linestyle_generator(
                colors=["black", "green", "red"],
                lines=['-',':'], markers=['o','s'], 
                hollow_styles=[False, False, True, True],
                )

    ax.plot(x, y, label="", 
            **linestyle.next())
    ax.set_xlim([0, 9])
    #ax.yaxis.set_ticks([0, 1, 2, 3, 4, 5, 6])

    ax.set_xlabel("z$_I$ ($\AA$)")
    ax.set_ylabel("Free Energy (kcal/mol)")
    ax.legend(loc='best')
    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)

plot_fig()


