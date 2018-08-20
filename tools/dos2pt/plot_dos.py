import os
import numpy as np
import peakutils
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpltex

PS2CM = 33.35641

def output_data(data, fname):
    o = open(fname, "w")
    for i in range(len(data[0])):
        for j in range(len(data)):
            o.write("%12.6f"%data[j][i])
        o.write("\n")
    o.close()

@mpltex.acs_decorator
def plot_dos():
    data = np.loadtxt("dos_sg")
    data = data.transpose()
    x = data[0] * PS2CM
    y1 = data[1]
    y2 = data[2]
    y3 = data[3]
    y4 = y1 + y2 + y3
    ymax = max(y4)
    y1 = y1/ymax
    y2 = y2/ymax
    y3 = y3/ymax

    peaks = [[],[]]
    #indexes = peakutils.indexes(y3, thres=0.02/max(y3), min_dist=100)
    indexes = peakutils.indexes(y3)
    for i in indexes:
        peaks[0].append(x[i])
        peaks[1].append(y3[i]+0.05)
        
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)

    #plt.scatter(peaks[0], peaks[1], marker="+", color="black")
    for i in range(len(indexes)):
        plt.text(peaks[0][i]-100, peaks[1][i], "%d"%peaks[0][i],)

    #get title
    title = os.getcwd().split("/")[-2]

    plt.xlim([-50,4000])
    plt.ylim([0,1])
    plt.xlabel("Wavenumber (cm-1)")
    plt.title(title)
    plt.tight_layout(pad=0.1)
    plt.savefig("dos-3-peaks.png", dpi=600)

@mpltex.acs_decorator
def plot_dos():
    data = np.loadtxt("dos_sg")
    data = data.transpose()
    x = data[0] * PS2CM
    y1 = data[1]
    y2 = data[2]
    y3 = data[3]
    y4 = y1 + y2 + y3

    plt.plot(x, y4)
    output_data([x, y4], "dos.dat")

    #get title
    title = os.getcwd().split("/")[-2]

    plt.xlim([-50,4000])
    #plt.ylim([0,1])
    plt.xlabel("Wavenumber (cm-1)")
    plt.title(title)
    plt.tight_layout(pad=0.1)
    plt.savefig("dos.png", dpi=600)


plot_dos()



