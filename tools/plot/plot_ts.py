# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_rxn(path):
    x0 = 0
    dx1 = path.dx1
    dx2 = path.dx2
    x1 = x0 + dx1[0]
    counter = 0
    for i in path.rxns:
        plt.plot([x0, x1], [i.e0, i.e0], lw=3, color=path.color)
        x0 = x1
        x1 = x0 + dx2[counter]
        if i.ea > 0:
            a = x0 + 0.5*dx2[counter]
            b = 4*i.ea/(dx2[counter]*dx2[counter])
            c = i.e0 + i.ea
            x = np.linspace(x0,a,9)
            y = -b*(x - a) * (x - a) + c          
            plt.plot(x,y, ls="--", color=path.color)
            a = x0 + 0.5*dx2[counter]
            b = 4*(i.ea + i.e0 - i.e1)/(dx2[counter]*dx2[counter])
            c = i.e0 + i.ea
            x = np.linspace(a,x1,9)
            y = -b*(x - a) * (x - a) + c          
            plt.plot(x,y, ls="--", color=path.color)
        else:
            plt.plot([x0, x1], [i.e0, i.e1], ls="--", color=path.color)
        x0 = x1
        x1 = x0 + dx1[counter+1]
        if counter == 0:
            plt.plot([x0, x1], [i.e1, i.e1], lw=3, 
                 color=path.color, label=path.name)
        else:
            plt.plot([x0, x1], [i.e1, i.e1], lw=3, 
                 color=path.color)
        counter += 1
    
class RxnPath():
    def __init__(self,):
        self.rxns = []
        self.color = "black"
        self.name = ""
        self.dx1 = [] # length for rxn and proc
        self.dx2 = [] # length for TS
        
class Rxn():
    def __init__(self,):
        self.e0 = 0.0
        self.e1 = 0.0
        self.ea = 0.0

def plot_path(pot, ts, dx1, dx2, color, name):
    for i in range(1, len(pot)):
        pot[i] = pot[i-1] + pot[i]
    path1 = RxnPath()
    for i in range(len(pot)-1):
        rxn = Rxn()
        rxn.e0 = pot[i]
        rxn.e1 = pot[i+1]
        rxn.ea = ts[i]    
        path1.rxns.append(rxn)
    path1.name = name
    path1.color = color
    path1.dx1 = dx1
    path1.dx2 = dx2
    
    plot_rxn(path1)

def main():
    pot = [0.000, -0.56, -0.48, -0.826, -1.478]
    ts = [0.0, 0.347, 0.068, 0.090]
    dx1 = [1.0]*len(pot)
    dx2 = [1.0]*len(ts)
    plot_path(pot, ts, dx1, dx2, "black", "reaxFF")
    
    pot = [0.0, -0.7, -0.70, -0.85, -2.0]
    ts = [0.0, 0.22, 0.0, 0.0]
    dx1 = [1.0]*len(pot)
    dx2 = [1.0]*len(ts)
    plot_path(pot, ts, dx1, dx2, "red", "QM")
    
    pot = [0.0, -0.247, 0.052, -0.048]
    ts = [0.088, 0.245, 0.283 ]
    dx1 = [1.0]*len(pot)
    dx2 = [3.0, 1.0, 1.0]
    plot_path(pot, ts, dx1, dx2, "blue", "reaxFF FE")
    
    plt.ylim([-5.0, 0.2])
    plt.ylabel("PE (ev)", size="x-large")
        
    frame = plt.gca()
    frame.axes.get_xaxis().set_ticks([])
    plt.legend(loc=3)
    plt.savefig("potential_map.png", dpi=600)
    plt.show()
    
if __name__ == "__main__":
    main()
