#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt

def parse_fort99():
    assert os.path.exists("fort.99")
    f = open("fort.99", 'r')
    reax = []
    qm = []

    for i in f:
        if i.strip().startswith("Energy"):
            reax.append(float(i[61:72]))
            qm.append(float(i[73:84]))
            break

    for i in f:
        if len(i.strip()) > 84:
            reax.append(float(i[61:72]))
            qm.append(float(i[73:84]))

    f.close()
    return reax, qm

def plot_all():
    xlabel = ''
    x = []
    if os.path.exists("bonds"):
        x = np.loadtxt("bonds")
        xlabel = r"bond length $\AA$"
    if os.path.exists("angles"):
        x = np.loadtxt("angles")
        xlabel = r"degree $^{\circ}$"
    reax, qm = parse_fort99()

    if len(x) == 0:
        x = np.linspace(0, 1, len(reax))

    assert len(x) == len(reax)

    plt.plot(x, reax, 'ro', ls='-', label="Reax")
    plt.plot(x, qm, 'bo', ls='-', label="QM")
    if xlabel == '':
        plt.xlabel(r"Unit-cell volume ($\AA ^3$)", size="x-large")
    else:
        plt.xlabel(xlabel, size="x-large")
    plt.ylabel(r"Potential energy (kcal)", size="x-large")
    plt.legend()
    plt.savefig("fort99.eps")
    plt.show()

def parse_trainset():
    f = open("trainset.in", 'r')
    for i in f:
        if i.strip().startswith("ENERGY"):
            break
    counter = 0
    title = ''
    names = []
    nframes = []
    for i in f:
        if i.strip().startswith("# "):
            if title:
                names.append(title)
                nframes.append(counter)
            title = i[1:].strip()
            counter = 0
        else:
            if i.strip().startswith("END"):
                pass
            else:
                counter +=1
    names.append(title)
    nframes.append(counter)
    f.close()

    return names, nframes

def plot_in_block():
    reax, qm = parse_fort99()
    names, nframes = parse_trainset()
    start = 0
    end = 0
    for i in range(len(names)):
        if i == 0:
            start = 0
            end = nframes[0]
        else:
            start += nframes[i-1]
            end += nframes[i]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(reax[start:end], 'ro', ls='-', label='ReaxFF')
        ax.plot(qm[start:end], 'bo', ls='--', label='QM')
        ax.set_title(names[i])
        ax.legend()
        plt.savefig("fig%02d.eps"%i)
        #latex
        tex = r"""\begin{figure}[htbp]
  \includegraphics[width=9cm]{fig%02d.eps}
  \caption{\label{fig:fig1}
  %s}
\end{figure}
"""%(i, names[i].replace('_', ' '))
        #print tex




if __name__ == "__main__":
    plot_all()
