import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def plot_out():
    data = np.loadtxt('out.dat')
    data = data.transpose()
    plt.plot(data[0], data[1])
    plt.show()
    print(data)

def plot_out_all():
    data = np.loadtxt('out1.dat')
    data = data.transpose()
    plt.errorbar(data[0], data[1]*0.043, data[2]*0.043, label='out1')

    data = np.loadtxt('out2.dat')
    data = data.transpose()
    plt.errorbar(data[0], data[1]*0.043, data[2]*0.043, label='out2')

    data = np.loadtxt('out3.dat')
    data = data.transpose()
    plt.errorbar(data[0], data[1]*0.043, data[2]*0.043, label='out3')
    plt.legend()
    plt.savefig('md.png')
    plt.show()
    print(data)

def plot_meta():
    fnames = []
    f = open('metadata', 'r')
    for i in f:
        tokens = i.strip().split()
        if len(tokens) == 3:
            fnames.append(tokens[0])
    f.close()
    
    for i in fnames:
        data = np.loadtxt(i)
        data = data.transpose()
        hist, bins = np.histogram(data[1], 20)
        print(hist, bins)
        mu = 100
        sigma = 15
        #y = mlab.normpdf(bins, mu, sigma)
        plt.plot(bins[1:], hist, '--')  
    plt.legend()
    plt.show()

#plot_meta()
plot_out_all()


