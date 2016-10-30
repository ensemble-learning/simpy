import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

@mpltex.acs_decorator
def plot_fig():
<<<<<<< HEAD
    data = np.loadtxt("all.xvg")
=======
    data = np.loadtxt("gofr.dat")
>>>>>>> 2fa3d11b619e5b1e1192c57a96f67798715b647c
    data = data.transpose()
    fig, ax = plt.subplots(1)
    linestyle = mpltex.linestyle_generator(
                colors=["black", "green", "red"],
                lines=['-',':'], markers=['o','s'], 
                hollow_styles=[False, False, True, True],
                )

<<<<<<< HEAD
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
=======
    x, y = [], []
    for i in range(len(data[0])):
        if i <= 50:
            x.append(data[0][i])
            y.append(data[1][i])
        else:
            if i % 10 == 0:
                x.append(data[0][i])
                y.append(data[1][i])
            
            
        
    ax.plot(x, y, 
            **linestyle.next())
    ax.set_xlabel("r$_{Li-O}$ ($\AA$)")
    ax.set_ylabel("g(r)")
    #ax.legend(loc='best')

    ax2 = ax.twinx()
    ax2.plot(data[0][::10], data[2][::10],
             **linestyle.next())
    ax2.plot([3.0, 8.0], [6.0, 6.0], ls="dotted", 
             lw=0.5, color="green")
    ax2.text(6.0, 6.2, "N = 6", color="green")
    ax2.set_ylim([0.0, 8.0])
    ax2.set_xlim([1.0, 8.0])
    ax2.set_ylabel("$\int$g(r) (N)", color="black")

>>>>>>> 2fa3d11b619e5b1e1192c57a96f67798715b647c
    fig.tight_layout(pad=0.1)
    fig.savefig("test.png", dpi=600)

plot_fig()
