import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_rdf_one(fname):
    label = fname.split('.')[0][1:]
    data = np.loadtxt(fname, skiprows=4)
    data = data.transpose()

    fig = plt.figure(figsize=(4, 3), dpi=300)
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    ax.plot(data[1], data[2], color='black')
    ax.scatter(data[1], data[2], 8, color='black')
    ax2.plot(data[1], data[3], '--', color='red')

    ax.set_xlim([1.0, 6.0])
    ax.set_xlabel('r($\AA$)', fontsize=14)
    ax.set_ylabel('g(r)', fontsize=14)
    ax.set_title(label, fontsize=14)

    ax2.set_ylabel('Coordination Number (N)', fontsize=14)
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    ax2.spines['right'].set_color('red')
    plt.tight_layout()
    print(fname.split('.'))
    plt.savefig(label+'.svg')
    plt.close()

if __name__ == '__main__':
    flist  = []
    f = open('flist', 'r')
    for i in f:
        flist.append(i.strip())
    f.close()
    for i in flist:
        plot_rdf_one(i)
