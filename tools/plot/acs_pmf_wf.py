import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
    data = np.loadtxt("all_wf.dat", 
                    dtype={'names': ('folder', 'x', 'dA/dx', 'err', 'A', 'wf'),
                         'formats': ('S20', 'f4', 'f4', 'f4', 'f4', 'f4')})
    xcut = 27.5 # free energy beyong 27.5 is too large
    x0 = 28.36  # the top layer of Cu
    
    x, y, wf = [], [], []
    for i in range(len(data)):
        tx = data[i][1]
        ty = data[i][4]
        twf = data[i][5]
        if tx < xcut:
            x.append(x0-tx)
            y.append(ty)
            wf.append(-twf)

    fig, ax = plt.subplots(1)
    linestyle = mpltex.linestyle_generator(
                colors=["black", "green", "red"],
                lines=['-',':'], markers=['o','s'], 
                hollow_styles=[False, False, True, True],
                )

    ax.plot(x, y, label="Free Energy", 
            **linestyle.next())
    #ax.yaxis.set_ticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_xlabel("z$_I$ ($\AA$)")
    ax.set_ylabel("Free Energy (eV)")
    ax.legend(loc='upper right', bbox_to_anchor=(0.88, 0.97))

    ax2 = ax.twinx()
    ax2.plot(x, wf, label="Work Function", 
            **linestyle.next())
    ax2.set_ylim([-5, -3])
    ax2.set_ylabel("Work Function (eV)", color="green")
    ax2.tick_params(axis='y', colors='green')
    ax2.legend(loc='upper right', bbox_to_anchor=(0.93, 0.87))

    ax.set_xlim([0, 9])
    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)

plot_fig()


