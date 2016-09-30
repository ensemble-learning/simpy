import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
    data = np.loadtxt("all.xvg")
    data = data.transpose()
    fig, ax = plt.subplots(1)
    linestyle = mpltex.linestyle_generator(
                colors=["black", "green", "red"],
                lines=['-',':'], markers=['o','s'], 
                hollow_styles=[False, False, True, True],
                )

    ax.plot(data[0][::10]*10, data[1][::10], label="k =  500", 
            **linestyle.next())
    ax.plot(data[2][::10]*10, data[3][::10] + 1.2, label="k = 1000", 
            **linestyle.next())
    ax.plot(data[4][::10]*10, data[5][::10] + 3.3, label="k = 1500", 
            **linestyle.next())
    ax.set_xlim([1.8, 6.7])
    ax.set_xlabel("r$_{Li-N}$ ($\AA$)")
    ax.set_ylabel("Free Energy (kcal/mol)")
    ax.legend(loc='best')
    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)

plot_fig()
